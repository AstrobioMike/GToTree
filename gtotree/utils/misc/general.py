import os
import sys
import gzip
import shutil
import json
import time
from dataclasses import dataclass, field, asdict, fields, is_dataclass
from typing import List
from tqdm import tqdm # type: ignore
import urllib.request
import urllib.error
from datetime import datetime
from gtotree.utils.misc.messaging import (report_message, color_text, wprint)
from gtotree.utils.misc.resume_state import invalidate_completed_from
from gtotree.utils.misc.stages import (GenomeRemovalStage,
                                  GENOME_REMOVAL_STAGE_ORDER,
                                  NCBI_REMOVAL_STAGES,
                                  PREPROCESSING_REMOVAL_STAGES,
                                  SCG_REMOVAL_STAGE_ORDER,
                                  PIPELINE_STAGES,
                                  stages_through,
                                  validate_stage)


def atomic_write_text(path, write_fn, encoding="utf-8"):
    tmp = f"{path}.part"
    try:
        with open(tmp, "w", encoding=encoding) as f:
            write_fn(f)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    except BaseException:
        try:
            os.remove(tmp)
        except OSError:
            pass
        raise


SEARCH_FAILURE_FILENAME = "SEARCH-FAILED.txt"


def record_search_failure(base_outpath, tool, detail):
    """
    Note why a per-genome search failed, in that genome's own results directory
    """
    try:
        os.makedirs(base_outpath, exist_ok=True)
        path = os.path.join(base_outpath, SEARCH_FAILURE_FILENAME)
        atomic_write_text(path, lambda f: f.write(f"{tool} failed: {detail}\n"))
    except BaseException:
        pass


def decode_pyhmmer_text(value):
    """
    Normalize a pyhmmer text attribute (`name`, `accession`, `description`) to str
    to be robust to a pyhmmer snag i hit earlier
    """
    if value is None:
        return None
    if isinstance(value, (bytes, bytearray)):
        return value.decode()
    return str(value)


@dataclass
class ToolsUsed:
    prodigal_used: bool = False
    gtdb_used: bool = False
    fasttree_used: bool = False
    veryfasttree_used: bool = False
    iqtree_used: bool = False
    pfam_db_used: bool = False
    kofamscan_used: bool = False
    universal_SCGs_used: bool = False


def download_with_tqdm(url, target, filename=None, leave=True,
                       retries=True, attempts=6, retry_wait=3,
                       speed_gate=False, min_mbps=2.0, probe_seconds=5.0,
                       probe_timeout=30):
    """
    Download `url` to `filename`, showing a tqdm progress bar, with optional
    transient-error retries and optional speed-gated route rerolling. Both are
    handled by a single unified attempt loop sharing one `attempts` budget.

    retries=True (default): transient failures (timeouts, connection resets,
    transient HTTP/URL errors) are retried up to `attempts` times, waiting
    `retry_wait` seconds between tries. A 404 is never retried (raised at once).
    Set retries=False for a single-shot download that raises on any error.

    speed_gate=False (default): don't judge throughput. speed_gate=True: GitHub
    release assets are served from a CDN whose per-connection throughput varies
    by edge node; a slow connection tends to stay slow, and reconnecting often
    lands on a faster edge. With the gate on, non-final attempts measure
    throughput over an initial `probe_seconds` window and, if it's below
    `min_mbps`, abort and reconnect to reroll the route. The final attempt
    accepts whatever speed it gets and runs to completion, so a persistently
    slow network still succeeds

    Raises the last underlying error if every attempt fails.
    """
    # a single attempt if retries are off
    if not retries:
        attempts = 1

    floor_bytes_per_s = (min_mbps * 1024 * 1024) if speed_gate else 0.0
    last_err = None

    for attempt in range(1, attempts + 1):
        is_final = (attempt == attempts)
        # enforce the speed floor only on non-final attempts (and only if gated)
        floor = 0.0 if is_final else floor_bytes_per_s
        desc = target

        try:
            _stream_once(url, filename, desc, leave, floor,
                         probe_seconds, connect_timeout=probe_timeout)
            if leave:
                sys.stderr.write("\n")
            return None

        except _TooSlow as e:
            # performance failure: reroll immediately, no wait
            last_err = e
            report_message(f"that was a slow route (try {attempt}/{attempts-1}), trying to get a faster one...", "yellow",
                           ii="          ", si="          ", width=90)
            print()
            continue

        except urllib.error.HTTPError as err:
            # a definitive 404 is never worth retrying
            if err.code == 404:
                raise
            last_err = err
            if is_final:
                raise
            wprint(color_text(
                f"    download failed (attempt {attempt}/{attempts}); retrying...",
                "yellow"))
            time.sleep(retry_wait)
            continue

        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as err:
            # transient network failure: wait, then retry
            last_err = err
            if is_final:
                raise
            wprint(color_text(
                f"    download failed (attempt {attempt}/{attempts}); retrying...",
                "yellow"))
            time.sleep(retry_wait)
            continue

    if last_err:
        raise last_err


class _TooSlow(Exception):
    def __init__(self, mbps):
        self.mbps = mbps
        super().__init__(f"too slow: {mbps:.2f} MB/s")


def _stream_once(url, filename, desc, leave, floor_bytes_per_s,
                 probe_seconds, connect_timeout=30):
    """
    Stream url->filename with a tqdm bar. The bar's total is resolved from the
    download response's own Content-Length header (no separate probe request),
    so a bounded bar shows whenever the server reports a length, and a slow
    route isn't paid for twice. If floor_bytes_per_s > 0, measure throughput
    over the first `probe_seconds` and raise _TooSlow if under floor.
    `connect_timeout` bounds the connect/read so a stalled server surfaces as a
    transient error rather than hanging. Network/HTTP errors propagate to the
    caller's attempt loop.
    """
    chunk = 1024 * 256  # 256 KB reads
    start = time.monotonic()
    probed = False
    req = urllib.request.Request(url, headers={'User-Agent': 'curl/8.0'})
    # Write to a temp file alongside the destination and atomically rename into
    # place only on a fully successful read. This way an interrupted download
    # (network drop, Ctrl-C, _TooSlow on the final attempt, process kill) never
    # leaves a truncated file at `filename` that a later os.path.isfile() check
    # would mistake for a complete download.
    tmp = f"{filename}.part"
    try:
        with urllib.request.urlopen(req, timeout=connect_timeout) as resp, open(tmp, "wb") as out, \
             tqdm(unit='B', unit_scale=True, unit_divisor=1024, miniters=1,
                  desc=desc, ncols=90, leave=leave, total=None) as bar:
            cl = resp.headers.get("Content-Length")
            total = int(cl) if cl else None
            if total is not None:
                bar.total = total
                bar.refresh()
            downloaded = 0
            while True:
                buf = resp.read(chunk)
                if not buf:
                    break
                out.write(buf)
                downloaded += len(buf)
                if total is not None:
                    bar.update(min(downloaded, total) - bar.n)
                else:
                    bar.update(len(buf))
                if floor_bytes_per_s > 0 and not probed:
                    elapsed = time.monotonic() - start
                    if elapsed >= probe_seconds:
                        probed = True
                        rate = downloaded / elapsed if elapsed > 0 else 0.0
                        if rate < floor_bytes_per_s:
                            bar.close()
                            raise _TooSlow(rate / (1024 * 1024))
            if total is not None and bar.n < total:
                bar.update(total - bar.n)
        os.replace(tmp, filename)
    except BaseException:
        # BaseException so Ctrl-C (KeyboardInterrupt) also cleans up the partial
        try:
            os.remove(tmp)
        except OSError:
            pass
        raise


GENOME_SOURCE_FIELDS = ("ncbi_accs", "genbank_files", "fasta_files", "amino_acid_files")

DERIVED_RUN_DATA_FIELDS = frozenset({"all_input_genomes"})


# How GenomeData.source values read in an output table. Shared so every table that
# reports where a genome came from says it the same way.
GENOME_SOURCE_LABELS = {
    "accession": "ncbi-accession",
    "genbank-file": "genbank-file",
    "nt-fasta-file": "nucleotide-fasta",
    "aa-fasta-file": "amino-acid-fasta",
}

# reference genomes pulled in by `-w` are accessions too, but they aren't accessions
# the USER listed, and a table that can't tell them apart can't answer "which of these
# did I ask for?" -- which is the whole point of the column
WANTED_REF_TAX_SOURCE_LABEL = "ncbi-accession (via -w)"


def genome_source_label(gd):
    """The output-table label for one genome's input source."""
    if getattr(gd, "from_wanted_ref_tax", False):
        return WANTED_REF_TAX_SOURCE_LABEL
    return GENOME_SOURCE_LABELS.get(gd.source, gd.source or "NA")


def genome_input_label(gd, run_data=None):
    """
    What the user actually handed over for this genome: the path they listed, or the
    `-w` request that produced it (which has no path to report).
    """
    if getattr(gd, "from_wanted_ref_tax", False):
        taxon = getattr(run_data, "wanted_ref_tax", None) if run_data else None
        return f"-w {taxon}" if taxon else "-w"
    return gd.provided_path or gd.id


@dataclass
class GenomeData:
    id: str
    source: str
    full_path: str
    provided_path: str
    basename: str
    taxid: str = None
    processing_done: bool = False
    final_AA_path: str = ""
    final_nt_path: str = ""
    prodigal_used: bool = False
    was_gzipped: bool = False
    acc_was_found: bool = None
    acc_was_downloaded: bool = None
    from_wanted_ref_tax: bool = False
    mapping: str = None
    hmm_search_done: bool = False
    hmm_search_failed: bool = None
    extract_seqs_failed: bool = None
    ko_search_done: bool = False
    ko_search_failed: bool = False
    pfam_search_done: bool = False
    pfam_search_failed: bool = False
    num_genes: int = None
    num_SCG_hits: int = None
    num_unique_SCG_hits: int = None
    num_SCG_hits_after_filtering: int = None
    reason_removed: str = None
    removed_at: str = None

    @property
    def removed(self):
        """
        Whether this genome is out of the run
        """
        return self.removed_at is not None

    @classmethod
    def from_path(cls, path, source):
        full_path = os.path.abspath(path)
        provided_path = path
        basename = os.path.basename(full_path)

        extensions_to_remove = [".gb", ".gbff", ".fasta", ".fna", ".fa", ".faa"]

        id = basename
        if id.lower().endswith(".gz"):
            id = id[:-3]

        for ext in extensions_to_remove:
            if id.lower().endswith(ext):
                id = id[:-len(ext)]
                break

        return cls(id=id,
                   source=source,
                   full_path=full_path,
                   provided_path=provided_path,
                   basename=basename)

    @classmethod
    def from_acc(cls, acc):
        full_path = None
        provided_path = None
        basename = acc
        id = acc
        source = "accession"

        return cls(id=id,
                   source=source,
                   full_path=full_path,
                   provided_path=provided_path,
                   basename=basename)

    def mark_processing_done(self, value=True):
        self.processing_done = value

    def mark_removed(self, reason, stage):
        """
        Drop this genome from the run

        `stage` is required and validated: a new removal site has to state where in the
        pipeline it sits, because every report about that stage is counted from it
        """
        validate_stage(stage, GENOME_REMOVAL_STAGE_ORDER, "genome-removal stage")
        self.removed_at = stage
        self.reason_removed = reason
        self.final_AA_path = None

    def mark_prodigal_used(self, value=True):
        self.prodigal_used = value

    def mark_was_gzipped(self, value=True):
        self.was_gzipped = value

    def mark_hmm_search_done(self, value=True):
        self.hmm_search_done = value
        self.hmm_search_failed = False

    def mark_hmm_search_failed(self, value=True):
        self.hmm_search_failed = value

    def mark_extract_seqs_failed(self, value=True):
        self.extract_seqs_failed = value

    def mark_ko_search_done(self, value=True):
        self.ko_search_done = value

    def mark_ko_search_failed(self, value=True):
        self.ko_search_done = True
        self.ko_search_failed = value

    def mark_pfam_search_done(self, value=True):
        self.pfam_search_done = value

    def mark_pfam_search_failed(self, value=True):
        self.pfam_search_done = True
        self.pfam_search_failed = value


@dataclass
class SCGset:
    id: str
    remaining: bool = None
    gene_length_filtered: bool = None
    num_genomes_with_any_hit: int = None            # incl. multi-copy hits
    num_genomes_with_hits: int = None               # contributing a usable seq
    num_genomes_with_hits_after_len_filtering: int = 0
    num_genomes_with_hits_after_genome_filtering: int = None
    aligned: bool = None
    trimmed: bool = None
    ready_for_cat: bool = None
    reason_removed: str = None
    removed_at: str = None

    @property
    def removed(self):
        return self.removed_at is not None

    @classmethod
    def from_id(cls, id):
        return cls(id, remaining=True, gene_length_filtered=False)

    def mark_removed(self, reason, stage):
        validate_stage(stage, SCG_REMOVAL_STAGE_ORDER, "SCG-removal stage")
        self.removed_at = stage
        self.reason_removed = reason
        self.remaining = False


@dataclass
class RunData:
    ncbi_accs: List[GenomeData] = field(default_factory=list)
    genbank_files: List[GenomeData] = field(default_factory=list)
    fasta_files: List[GenomeData] = field(default_factory=list)
    amino_acid_files: List[GenomeData] = field(default_factory=list)
    all_input_genomes: List[GenomeData] = field(default_factory=list)
    SCG_targets: List[SCGset] = field(default_factory=list)

    start_time: datetime = None
    ncbi_sub_table_path: str = ""
    ncbi_processing_dir: str = ""
    ncbi_processing_dir_rel: str = ""
    genbank_processing_dir: str = ""
    fasta_processing_dir: str = ""
    AA_processing_dir: str = ""
    output_dir: str = ""
    output_dir_rel: str = ""
    run_files_dir: str = ""
    run_files_dir_rel: str = ""
    individual_gene_alignments_dir: str = ""
    individual_gene_alignments_dir_rel: str = ""
    run_data_path: str = ""
    hmm_path: str = ""
    mapping_file_path: str = ""
    mapping_dict: dict = field(default_factory=dict)
    tax_info_dict: dict = field(default_factory=dict)
    initial_mapping_IDs_from_user: List[str] = field(default_factory=list)
    ready_genome_files_dir: str = ""
    hmm_results_dir: str = ""
    found_SCG_seqs_dir: str = ""
    tmp_dir: str = ""
    log_file: str = ""
    logs_dir: str = ""
    logs_dir_rel: str = ""
    gtotree_logs_dir: str = ""
    best_hit_mode: bool = False
    seq_length_cutoff: float = None
    gene_representation_cutoff: float = None
    completed_stages: dict = field(default_factory=dict)
    updating_headers: bool = False
    use_muscle_super5: bool = False
    num_muscle_threads: int = 5
    nucleotide_mode: bool = False
    general_ext: str = ""
    concatenated_alignment_path: str = ""
    final_alignment_path: str = ""
    final_alignment_length: int = 0
    original_tree_path: str = ""
    final_tree_path: str = ""
    target_kos_file: str = None
    total_ko_targets: int = 0
    target_kos_tsv: str = None
    target_ko_profiles_dir: str = None
    target_pfams_file: str = None
    total_pfam_targets: int = 0
    all_pfam_targets_hmm_path: str = None
    additional_pfam_searching_done: bool = False
    additional_ko_searching_done: bool = False

    pfam_dict: dict = field(default_factory=dict)
    pfam_results_dir: str = ""
    pfam_results_dir_rel: str = ""
    tmp_pfam_results_dir: str = ""
    wanted_pfam_targets: List[str] = field(default_factory=list)
    found_pfam_targets: List[str] = field(default_factory=list)
    failed_pfam_targets: List[str] = field(default_factory=list)

    ko_results_dir: str = ""
    ko_results_dir_rel: str = ""
    tmp_ko_results_dir: str = ""
    wanted_ko_targets: List[str] = field(default_factory=list)
    found_ko_targets: List[str] = field(default_factory=list)
    failed_ko_targets: List[str] = field(default_factory=list)

    tools_used: ToolsUsed = field(default_factory=ToolsUsed)

    # --- how a `-w` selection was made ---------------------------------------
    # Kept on run_data rather than passed around as the RefGenomeSelection object so
    # the reporting can describe the selection on a resumed run too, where `-w` isn't
    # re-resolved and there is no selection object to ask. Plain str/int so they
    # survive the run-data.json round-trip.
    wanted_ref_tax: str = None
    wanted_ref_tax_rank: str = None
    wanted_ref_tax_derep_rank: str = None
    wanted_ref_tax_num_selected: int = 0

    # dict of the run parameters that affect what this run produces; compared on -R
    # to decide whether resuming is safe. See preflight_checks.build_fingerprint
    fingerprint: dict = field(default_factory=dict)

    # --- pipeline-stage completion --------------------------------------------

    def mark_stage_complete(self, stage):
        validate_stage(stage, PIPELINE_STAGES, "pipeline stage")
        self.completed_stages[stage] = datetime.now().isoformat(timespec="seconds")

    def stage_is_complete(self, stage) -> bool:
        validate_stage(stage, PIPELINE_STAGES, "pipeline stage")
        return stage in self.completed_stages

    def invalidate_stages_from(self, stage):
        """
        Forget `stage` and everything downstream of it

        Anything computed from a stage's output stops being trustworthy the moment that
        stage has to be redone, even if it was marked complete
        """
        validate_stage(stage, PIPELINE_STAGES, "pipeline stage")
        invalidate_completed_from(self.completed_stages, stage, PIPELINE_STAGES)

    # --- stage-scoped removals ------------------------------------------------

    def genomes_removed_at(self, *stages, source=None) -> List[GenomeData]:
        """
        Input genomes dropped at any of `stages`

        `source` optionally narrows to one of GENOME_SOURCE_FIELDS, so a per-source
        report doesn't count another source's failures
        """
        for stage in stages:
            validate_stage(stage, GENOME_REMOVAL_STAGE_ORDER, "genome-removal stage")
        wanted = frozenset(stages)
        pool = getattr(self, source) if source else self.all_input_genomes
        return [gd for gd in pool if gd.removed_at in wanted]

    def genomes_alive_through(self, stage, source=None) -> List[GenomeData]:
        """
        Input genomes that were still in the run as of the end of `stage`

        This is the count a stage's report wants. `get_all_remaining_input_genomes()`
        answers the different question "who's left *now*?" (which wasn't helpful for resuming counts)
        """
        by_now = stages_through(stage, GENOME_REMOVAL_STAGE_ORDER)
        pool = getattr(self, source) if source else self.all_input_genomes
        return [gd for gd in pool if gd.removed_at not in by_now]

    def SCG_targets_removed_at(self, *stages) -> List[SCGset]:
        for stage in stages:
            validate_stage(stage, SCG_REMOVAL_STAGE_ORDER, "SCG-removal stage")
        wanted = frozenset(stages)
        return [scg for scg in self.SCG_targets if scg.removed_at in wanted]

    def SCG_targets_alive_through(self, stage) -> List[SCGset]:
        by_now = stages_through(stage, SCG_REMOVAL_STAGE_ORDER)
        return [scg for scg in self.SCG_targets if scg.removed_at not in by_now]

    def get_all_SCG_targets(self) -> List[SCGset]:
        return list(self.SCG_targets)

    def get_all_SCG_targets_remaining(self) -> List[SCGset]:
        return [scg for scg in self.SCG_targets if scg.remaining and not scg.removed]

    def get_all_SCG_targets_remaining_but_not_filtered(self) -> List[SCGset]:
        return [scg for scg in self.SCG_targets if scg.remaining and not scg.gene_length_filtered]

    def get_all_SCG_targets_ready_for_concatenation(self) -> List[SCGset]:
        return [scg for scg in self.SCG_targets if scg.ready_for_cat and not scg.removed]

    def update_all_input_genomes(self):
        """
        Build `all_input_genomes` from the four source lists
        """
        self.all_input_genomes = []
        for name in GENOME_SOURCE_FIELDS:
            self.all_input_genomes.extend(getattr(self, name))

    def get_all_input_genome_ids(self) -> List[str]:
        return [gd.id for gd in self.all_input_genomes]

    def merge_wanted_ref_tax_accessions(self, accessions) -> int:
        """
        Fold --wanted-ref-tax (-w) accessions into the NCBI-accession input pool,
        skipping any already provided by the user by accession id.
        Order is preserved: existing accs first, then the new -w ones in the order
        the taxonomy core returned them. Refreshes all_input_genomes and returns the
        number actually added.
        """
        existing = {gd.id for gd in self.ncbi_accs}
        added = 0
        for acc in accessions:
            if acc in existing:
                continue
            gd = GenomeData.from_acc(acc)
            gd.from_wanted_ref_tax = True
            self.ncbi_accs.append(gd)
            existing.add(acc)
            added += 1
        if added:
            self.update_all_input_genomes()
        return added

    def record_wanted_ref_tax_selection(self, selection, taxon=None):
        """
        Remember HOW a `-w` selection was made, for the reporting layer.

        `wanted_ref_tax_num_selected` is what the taxonomy core handed back, BEFORE
        deduping against `-a`, so the difference between it and
        `get_wanted_ref_tax_accs()` is exactly the overlap worth mentioning.
        """
        if selection is None:
            return self

        self.wanted_ref_tax = taxon or getattr(selection, "canonical", None)
        self.wanted_ref_tax_rank = selection.resolved_rank
        self.wanted_ref_tax_derep_rank = selection.effective_derep_rank
        self.wanted_ref_tax_num_selected = len(selection.accessions)

        return self

    def get_all_processed_genomes(self) -> List[GenomeData]:
        return [gd for gd in self.all_input_genomes if gd.processing_done]

    def get_all_input_genomes_for_filtering(self) -> List[GenomeData]:
        return [gd for gd in self.all_input_genomes if gd.processing_done and gd.hmm_search_done and not gd.removed]

    def get_all_input_genomes_due_for_SCG_min_hit_filtering(self) -> List[GenomeData]:
        # was a match against two literal `reason_removed` strings, which would have
        # silently started returning [] the moment either bit of wording was edited
        return self.genomes_removed_at(GenomeRemovalStage.SCG_HIT_FILTER)

    def get_input_ncbi_accs(self) -> List[str]:
        return [gd.id for gd in self.ncbi_accs]

    def get_user_provided_ncbi_accs(self) -> List[str]:
        return [gd.id for gd in self.ncbi_accs if not gd.from_wanted_ref_tax]

    def get_wanted_ref_tax_accs(self) -> List[str]:
        return [gd.id for gd in self.ncbi_accs if gd.from_wanted_ref_tax]

    def get_input_genbank_ids(self) -> List[str]:
        return [gd.id for gd in self.genbank_files]

    def get_input_fasta_ids(self) -> List[str]:
        return [gd.id for gd in self.fasta_files]

    def get_input_amino_acid_ids(self) -> List[str]:
        return [gd.id for gd in self.amino_acid_files]

    def remaining_ncbi_accs(self) -> List[str]:
        return [gd.id for gd in self.ncbi_accs if not gd.removed]

    def get_ncbi_accs_not_downloaded(self) -> List[str]:
        return [gd.id for gd in self.ncbi_accs if gd.acc_was_downloaded is False]

    def get_ncbi_accs_not_found(self) -> List[str]:
        return [gd.id for gd in self.ncbi_accs if gd.acc_was_found is False]

    def get_removed_ncbi_accs(self) -> List[str]:
        """Accessions dropped *during NCBI handling*"""
        return [gd.id for gd in
                self.genomes_removed_at(*NCBI_REMOVAL_STAGES, source="ncbi_accs")]

    def get_genomes_removed_during_processing(self) -> List[GenomeData]:
        return self.genomes_removed_at(*PREPROCESSING_REMOVAL_STAGES)

    def get_all_removed_input_genomes(self) -> List[str]:
        return [gd.id for gd in self.all_input_genomes if gd.removed]

    def get_all_remaining_input_genomes(self) -> List[GenomeData]:
        return [gd for gd in self.all_input_genomes if not gd.removed]

    def get_all_remaining_input_genome_ids(self) -> List[str]:
        return [gd.id for gd in self.all_input_genomes if not gd.removed]

    def get_done_ncbi_accs(self) -> List[GenomeData]:
        return [gd for gd in self.ncbi_accs if gd.processing_done]

    def _failed_at(self, stage, source) -> List[GenomeData]:
        return self.genomes_removed_at(stage, source=source)

    def get_failed_genbank_ids(self) -> List[str]:
        return [gd.id for gd in
                self._failed_at(GenomeRemovalStage.GENBANK_PREP, "genbank_files")]

    def get_failed_fasta_ids(self) -> List[str]:
        return [gd.id for gd in
                self._failed_at(GenomeRemovalStage.FASTA_PREP, "fasta_files")]

    def get_failed_amino_acid_ids(self) -> List[str]:
        return [gd.id for gd in
                self._failed_at(GenomeRemovalStage.AMINO_ACID_PREP, "amino_acid_files")]

    def get_prodigal_used_genbank_ids(self) -> List[str]:
        return [gd.id for gd in self.genbank_files if gd.prodigal_used]


def populate_run_data(args):
    run_data = RunData()

    if args.ncbi_accessions:
        with open(args.ncbi_accessions) as f:
            entries_list = f.read().splitlines()
        run_data.ncbi_accs = [GenomeData.from_acc(entry) for entry in entries_list]

    if args.genbank_files:
        with open(args.genbank_files) as f:
            entries_list = f.read().splitlines()
        run_data.genbank_files = [GenomeData.from_path(entry, "genbank-file") for entry in entries_list]

    if args.fasta_files:
        with open(args.fasta_files) as f:
            entries_list = f.read().splitlines()
        run_data.fasta_files = [GenomeData.from_path(entry, "nt-fasta-file") for entry in entries_list]
        for gd in run_data.fasta_files:
            gd.prodigal_used = True

    if args.amino_acid_files:
        with open(args.amino_acid_files) as f:
            entries_list = f.read().splitlines()
        run_data.amino_acid_files = [GenomeData.from_path(entry, "aa-fasta-file") for entry in entries_list]

    run_data.update_all_input_genomes()
    run_data.run_files_dir = args.run_files_dir
    run_data.run_files_dir_rel = args.run_files_dir_rel
    run_data.run_data_path = run_data.run_files_dir + "/run-data.json"

    return run_data


def run_data_as_dict(run_data):
    """
    Build the JSON-serializable view of a RunData, skipping DERIVED_RUN_DATA_FIELDS
    """
    run_data_dict = {}
    for f in fields(run_data):
        if f.name in DERIVED_RUN_DATA_FIELDS:
            continue
        value = getattr(run_data, f.name)
        if is_dataclass(value):
            value = asdict(value)
        elif isinstance(value, list):
            value = [asdict(item) if is_dataclass(item) else item for item in value]
        run_data_dict[f.name] = value
    return run_data_dict


def write_run_data(run_data):
    run_data_dict = run_data_as_dict(run_data)
    if isinstance(run_data.start_time, datetime):
        run_data_dict['start_time'] = run_data.start_time.isoformat()
    atomic_write_text(run_data.run_data_path,
                      lambda f: json.dump(run_data_dict, f, indent=2))


class CorruptRunData(Exception):
    """The saved run-data.json exists but can't be parsed back into a RunData."""


def _known_fields_only(cls, values):
    """
    Drop keys that aren't fields of `cls`

    A run-data.json written by an older GToTree can carry fields that no longer exist
    (`processing_failed`, `run_complete`, ...). Without this, loading one raises
    TypeError and surfaces as "corrupt run data". Letting it instead load means the
    fingerprint's `state_version` gets to do its job and refuse the resume with an
    explanation instead.
    """
    allowed = {f.name for f in fields(cls)}
    return {k: v for k, v in values.items() if k in allowed}


def read_run_data(path):
    """
    Load a saved RunData, or return None if there isn't one at `path`
    """
    try:
        with open(path) as f:
            run_data_dict = json.load(f)

        for name in DERIVED_RUN_DATA_FIELDS:
            run_data_dict.pop(name, None)

        for name in GENOME_SOURCE_FIELDS:
            if name in run_data_dict:
                run_data_dict[name] = [
                    GenomeData(**_known_fields_only(GenomeData, gd))
                    if isinstance(gd, dict) else gd
                    for gd in run_data_dict[name]]

        if "tools_used" in run_data_dict and run_data_dict["tools_used"] is not None:
            run_data_dict["tools_used"] = ToolsUsed(
                **_known_fields_only(ToolsUsed, run_data_dict["tools_used"]))
        else:
            run_data_dict["tools_used"] = ToolsUsed()

        if isinstance(run_data_dict.get("start_time"), str):
            run_data_dict["start_time"] = datetime.fromisoformat(run_data_dict["start_time"])

        if "SCG_targets" in run_data_dict:
            run_data_dict["SCG_targets"] = [
                SCGset(**_known_fields_only(SCGset, scg))
                if isinstance(scg, dict) else scg
                for scg in run_data_dict["SCG_targets"]]

        run_data = RunData(**_known_fields_only(RunData, run_data_dict))
        run_data.update_all_input_genomes()

        return run_data

    except FileNotFoundError:
        return None
    except (json.JSONDecodeError, TypeError, ValueError) as e:
        raise CorruptRunData(
            f"the saved run data at \"{path}\" couldn't be read back ({e})") from e


def touch(path):
    with open(path, 'a'):
        os.utime(path, None)


class OutputDirExistsError(Exception):
    """
    The output directory already exists and neither `-R` nor `-F` was given
    """


def prepare_output_dir(out_dir, resume, force_overwrite,
                       work_dir_name="working-dir", ii="  ", si="  "):
    """
    Create the output dir and its working dir, honoring `-F` and `-R`
    """
    out_dir = out_dir.rstrip("/")
    work_dir = os.path.join(out_dir, work_dir_name)

    if os.path.exists(out_dir):
        if resume:
            if not os.path.isdir(work_dir):
                report_message(
                    f"There's no working directory in '{out_dir}' from a previous run "
                    "to resume from, so we'll start fresh.", "yellow", ii=ii, si=si)
            os.makedirs(work_dir, exist_ok=True)
            return out_dir, work_dir

        if not force_overwrite:
            raise OutputDirExistsError(
                f"The output directory '{out_dir}' already exists, and we don't want "
                "to overwrite anything accidentally. Use `-R` to resume that run, "
                "`-F` to overwrite it, or specify a different directory with `-o`.")
        shutil.rmtree(out_dir)

    elif resume:
        report_message(
            f"`-R`/`--resume` was specified, but '{out_dir}' doesn't exist yet, so "
            "we'll just start fresh.", "yellow", ii=ii, si=si)

    os.makedirs(work_dir, exist_ok=True)

    return out_dir, work_dir


GTT_PROGRESS_BAR_FORMAT = (
    "      {percentage:3.0f}%|{bar}| "
    "{n_fmt}/{total_fmt} "
    "[time elapsed: {elapsed} | est. remaining: {remaining}]"
)

GTT_PROGRESS_BAR_FORMAT_NO_COUNT_INDENTED = (
    "        {percentage:3.0f}%|{bar}| "
    "[time elapsed: {elapsed} | est. remaining: {remaining}]"
)

GTT_PROGRESS_BAR_FORMAT_INDENTED = (
    "        {percentage:3.0f}%|{bar}| "
    "{n_fmt}/{total_fmt} "
    "[time elapsed: {elapsed} | est. remaining: {remaining}]"
)

GTT_PROGRESS_BAR_FORMAT_INDENTED_6 = (
    "      {percentage:3.0f}%|{bar}| "
    "{n_fmt}/{total_fmt} "
    "[time elapsed: {elapsed} | est. remaining: {remaining}]"
)


# `smoothing=0` makes tqdm use a plain cumulative average (n / elapsed) instead, which is
# what i want in almost everything
GTT_PROGRESS_SMOOTHING = 0


# How many items to keep queued per worker thread in run_pooled_stage. Enough that a
# thread finishing early always has work waiting, small enough that the queued futures
# and their results don't dominate memory on a large run
WORKER_QUEUE_DEPTH = 4


def required_count(total, fraction):
    """
    How many of `total` are needed to meet `fraction`, rounding halves up
    """
    return int(total * fraction + 0.5)


def search_threads_per_genome(args):
    """
    Threads to hand an individual per-genome search (hmmsearch/kofamscan)

    i'm keeping this as a function for now because i keep changing my mind if i
    want it to be constant to determined...
    """
    return 2


def run_pooled_stage(items, worker, apply_result, args, run_data,
                     max_workers_cap=None, bar_format=None, lead_newline=True,
                     smoothing=GTT_PROGRESS_SMOOTHING):
    """
    Dispatches `worker(item, run_data)` across a ThreadPoolExecutor and, as each
    future completes, calls `apply_result(item, result, run_data)` on the MAIN
    thread. This split is deliberate and load-bearing:

      - `worker` runs concurrently in a thread. It must only touch per-item files
        and return a plain result object; it must NOT mutate shared run_data
        state or append to shared output files (that would race)
      - `apply_result` runs single-threaded here, in completion order, so all
        GenomeData/SCG bookkeeping, shared-table appends, and shared counters are
        safe

    `max_workers_cap` optionally clamps the worker count below --num-jobs (e.g., like
    how i limit concurrent downloads from ncbi); None means just use --num-jobs

    `bar_format` overrides the tqdm bar format (default: GTT_PROGRESS_BAR_FORMAT).

    `lead_newline` prints a blank line before the bar (the default). I needed to add this
    after i started formatting differently in some uses

    `smoothing` is handed straight to tqdm and defaults to GTT_PROGRESS_SMOOTHING (0,
    i.e., a whole-run average), so it's more informative in this situation

    NOTE ON WORKERS AND EXCEPTIONS: `worker` must not raise. A worker exception
    propagates out of `future.result()` here and aborts the whole stage partway
    through, with some items applied and some not. Workers should wrap their body in
    `try/except BaseException` and return a status object describing the failure, so
    `apply_result` can record it and the stage can carry on. All the processing and
    search workers follow this.

    On KeyboardInterrupt, queued-but-not-yet-started work is cancelled rather than
    drained, so a ctrl-c during a big stage exits once the in-flight workers finish
    instead of waiting for every remaining item.
    """
    from concurrent.futures import ThreadPoolExecutor, FIRST_COMPLETED, wait
    from contextvars import copy_context
    from itertools import islice

    num_workers = max(1, args.num_jobs)
    if max_workers_cap is not None:
        num_workers = min(num_workers, max_workers_cap)

    # Only ever keep this many items submitted at once, rather than queueing all of
    # them up front (saves storage and time on an interrupt)
    window = max(num_workers * WORKER_QUEUE_DEPTH, num_workers + 1)

    if lead_newline:
        print("")
    pool = ThreadPoolExecutor(max_workers=num_workers)

    def submit(item):
        """
        Submit `worker` carrying a copy of this thread's context.

        A new thread starts with an empty context, so ContextVars set on the main
        thread fall back to their defaults inside a worker. That matters for
        `log_file_var`: preflight sets it to an absolute path inside the output
        directory, but a worker would see the relative default and any reporter it
        called would append `gtotree-runlog.txt` to the user's cwd instead. No worker
        currently reaches a reporter, so this defends the invariant rather than fixing
        a live bug -- but it's an easy one to break by adding a per-item warning.

        A fresh copy per submission: a Context cannot be entered re-entrantly, so
        reusing one across concurrent workers would raise.
        """
        return pool.submit(copy_context().run, worker, item, run_data)

    remaining = iter(items)
    in_flight = {}
    try:
        with tqdm(total=len(items),
                  bar_format=bar_format or GTT_PROGRESS_BAR_FORMAT, ncols=76,
                  smoothing=smoothing) as pbar:

            for item in islice(remaining, window):
                in_flight[submit(item)] = item

            while in_flight:
                done, _ = wait(in_flight, return_when=FIRST_COMPLETED)

                for future in done:
                    item = in_flight.pop(future)
                    result = future.result()
                    # apply_result still runs here on the MAIN thread, in completion
                    # order, exactly as before
                    apply_result(item, result, run_data)
                    pbar.update(1)

                # top up by however many just completed, so the window stays full
                for item in islice(remaining, len(done)):
                    in_flight[submit(item)] = item

    except KeyboardInterrupt:
        pool.shutdown(wait=False, cancel_futures=True)
        raise
    finally:
        pool.shutdown(wait=True)

    return run_data


def gunzip_if_needed(path):
    if path.endswith(".gz"):
        gunzipped_path = path[:-3]
        with gzip.open(path, 'rb') as f_in, open(gunzipped_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        return gunzipped_path, True
    return path, False


def remove_file_if_exists(path):
    try:
        os.remove(path)
    except FileNotFoundError:
        pass


def file_is_usable_else_clear(path):
    """
    True if `path` exists and is non-empty. **Deletes it** if it exists but is empty
    """
    try:
        if os.path.getsize(path) > 0:
            return True
        remove_file_if_exists(path)
    except FileNotFoundError:
        pass
    return False


def concat_files(file_list, output_file):
    """
    Concatenate `file_list` into `output_file`, atomically
    """
    def _write(outfile):
        for fname in file_list:
            with open(fname) as infile:
                shutil.copyfileobj(infile, outfile)

    atomic_write_text(output_file, _write)


def cleanup(args, run_data):
    if not args.debug:
        shutil.rmtree(run_data.tmp_dir, ignore_errors=True)
        run_data.tmp_dir = ""
