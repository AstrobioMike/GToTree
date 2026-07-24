"""
Local genome-file input for `gtt gen-scg-hmms`

The main GToTree driver accepts four genome sources: NCBI accessions (`-a`), GenBank
files (`-g`), nucleotide fastas (`-f`), and amino-acid fastas (`-A`). This module brings
the three local-file types to gen-scg-hmms so it can accept the same inputs.

Rather than reimplement the parsing, I'm re-using the main run's actual preprocessing
logic from `gtotree.utils.seqs`:

    extract_filter_and_rename_cds_amino_acids_from_gb   GenBank CDS translations
    extract_fasta_from_gb                               GenBank -> nucleotide fasta
    _filter_and_rename_fasta                            length filter + header renaming

Those functions take a `run_data` but only ever read a few directory attributes off it
(`ready_genome_files_dir`, `genbank_processing_dir`, `nucleotide_mode`).
So `_ShimRunData` supplies what's needed, letting this standalone helper
share the real logic instead of reimplementing more

Genome IDs come from `GenomeData.from_path()`, the same factory the main driver uses, so
a given file yields the same ID in both programs (extension and `.gz` stripping included).
"""

import os
import gzip
import shutil

from gtotree.utils.general import GenomeData, remove_file_if_exists
from gtotree.utils.seqs import (extract_filter_and_rename_cds_amino_acids_from_gb,
                                extract_fasta_from_gb,
                                _filter_and_rename_fasta)
from gtotree.utils.hmms.gen_scg_hmms import GenSCGHMMsError
from gtotree.utils.hmms.gen_scg_hmms_genomes import (TargetGenomeError,
                                                     MISSED_NO_PROTEINS,
                                                     MISSED_PRODIGAL_FAILED,
                                                     run_prodigal)


# sources, matching the main driver's naming
SOURCE_GENBANK = "genbank"
SOURCE_FASTA = "fasta"
SOURCE_AMINO_ACID = "amino-acid"

MISSED_UNREADABLE = "file could not be read"
MISSED_NO_CDS = "no usable CDS translations found"


class _ShimRunData:
    """
    Minimal stand-in for RunData for the `seqs.py` helpers.

    Those helpers only read directory paths and `nucleotide_mode` off run_data, so this
    supplies just those
    """

    __slots__ = ("ready_genome_files_dir", "genbank_processing_dir",
                 "fasta_processing_dir", "ncbi_processing_dir", "nucleotide_mode")

    def __init__(self, work_dir):
        self.ready_genome_files_dir = work_dir
        self.genbank_processing_dir = work_dir
        self.fasta_processing_dir = work_dir
        self.ncbi_processing_dir = work_dir
        self.nucleotide_mode = False


def read_paths_file(path, label):
    """
    Read a file listing genome-file paths, one per line.

    Matches the main driver's input convention: `-g`, `-f`, and `-A` each take a file
    *listing* paths, not the genome files themselves. Blank lines and `#` comments are
    ignored, order is preserved, and duplicates are dropped.
    """
    if not os.path.isfile(path):
        raise GenSCGHMMsError(
            f"the specified {label} file '{path}' can't be found.")

    paths = []
    seen = set()
    with open(path) as f:
        for line in f:
            entry = line.strip()
            if not entry or entry.startswith("#"):
                continue
            if entry not in seen:
                seen.add(entry)
                paths.append(entry)

    if not paths:
        raise GenSCGHMMsError(f"no paths were read out of '{path}'.")

    return paths


def build_local_genomes(args):
    """
    Turn the -g/-f/-A inputs into GenomeData objects.

    Returns (genomes, missing) where `missing` holds (id, reason) for any listed path
    that doesn't exist on disk, checked up front so a typo in a long list is reported
    before any expensive work starts, rather than midway through.
    """
    genomes = []
    missing = []

    for attr, label, source in (
        ("genbank_files", "genbank-files", SOURCE_GENBANK),
        ("fasta_files", "fasta-files", SOURCE_FASTA),
        ("amino_acid_files", "amino-acid-files", SOURCE_AMINO_ACID),
    ):
        listing = getattr(args, attr, None)
        if not listing:
            continue

        for path in read_paths_file(listing, label):
            gd = GenomeData.from_path(path, source)
            if not os.path.isfile(gd.full_path):
                missing.append((gd.id, f"{source} file not found: {path}"))
                continue
            genomes.append(gd)

    return genomes, missing


def _gunzip_into_workdir(gd, work_dir):
    """
    Decompress a gzipped input into the working directory.
    """
    if not gd.full_path.endswith(".gz"):
        return gd.full_path, False

    inner_name = os.path.basename(gd.full_path)[:-3]
    staged = os.path.join(work_dir, f"_unzipped_{gd.id}_{inner_name}")

    with gzip.open(gd.full_path, "rb") as f_in, open(staged, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out, length=1024 * 1024)

    return staged, True


def process_local_genome(gd, work_dir, nucleotide_fallback=True):
    """
    Produce an amino-acid fasta for one local genome file

    Dispatches on `gd.source`:

      amino-acid  already proteins -- just length-filter and rename headers
      genbank     try CDS translations first, fall back to nucleotides + prodigal
      fasta       nucleotides, so always prodigal

    Returns (aa_path, used_prodigal). Raises TargetGenomeError with a report-ready reason
    on failure. Gzipped inputs are staged into the working dir and cleaned up; the user's
    input files are never modified or removed.
    """
    shim = _ShimRunData(work_dir)
    path = gd.full_path
    needs_cleanup = False

    try:
        path, needs_cleanup = _gunzip_into_workdir(gd, work_dir)

        if gd.source == SOURCE_AMINO_ACID:
            return _process_amino_acids(gd, path, shim), False

        if gd.source == SOURCE_GENBANK:
            return _process_genbank(gd, path, shim, work_dir, nucleotide_fallback)

        if gd.source == SOURCE_FASTA:
            return _process_fasta(gd, path, work_dir), True

        raise TargetGenomeError(f"unrecognized genome source '{gd.source}'")

    except TargetGenomeError:
        raise
    except OSError as e:
        raise TargetGenomeError(f"{MISSED_UNREADABLE}: {e}")
    except Exception as e:  # noqa: BLE001 - surfaced per-genome, run continues
        raise TargetGenomeError(f"unexpected failure: {e}")
    finally:
        if needs_cleanup:
            remove_file_if_exists(path)


def _process_amino_acids(gd, path, shim):
    """Length-filter and rename an existing amino-acid fasta."""
    out_path = os.path.join(shim.ready_genome_files_dir, f"{gd.id}.faa")
    num = _filter_and_rename_fasta(out_path, path, gd.id, 99999)

    if num == 0:
        remove_file_if_exists(out_path)
        raise TargetGenomeError(MISSED_NO_PROTEINS)

    return out_path


def _process_genbank(gd, path, shim, work_dir, nucleotide_fallback):
    """
    GenBank: prefer the CDS translations already in the file if present
    Fall back to prodigal when there are none usable
    """
    done, aa_path, num = extract_filter_and_rename_cds_amino_acids_from_gb(
        gd.id, path, shim)

    if done and num > 0:
        return aa_path, False

    if not nucleotide_fallback:
        raise TargetGenomeError(MISSED_NO_CDS)

    # no usable translations -> pull the nucleotides out and call genes
    extract_fasta_from_gb(gd.id, path, shim)
    nt_path = os.path.join(shim.genbank_processing_dir, f"{gd.id}.fasta")

    if not (os.path.isfile(nt_path) and os.path.getsize(nt_path) > 0):
        raise TargetGenomeError(MISSED_NO_CDS)

    try:
        aa_path = _prodigal_to_ready(gd, nt_path, shim, work_dir)
    finally:
        remove_file_if_exists(nt_path)

    return aa_path, True


def _process_fasta(gd, path, work_dir):
    """Nucleotide fasta: always needs gene calling."""
    shim = _ShimRunData(work_dir)
    return _prodigal_to_ready(gd, path, shim, work_dir)


def _prodigal_to_ready(gd, nt_path, shim, work_dir):
    """
    Call genes with prodigal, then length-filter and rename into the ready file.

    Reuses this helper's own `run_prodigal` (atomic, raises on failure) rather than the
    main gtotree program's (because that one writes into stage-specific directories and
    reports through run_data
    """
    raw_aa_path = os.path.join(work_dir, f"{gd.id}_prodigal.faa")

    try:
        run_prodigal(nt_path, raw_aa_path)
    except TargetGenomeError:
        raise
    except Exception as e:
        raise TargetGenomeError(f"{MISSED_PRODIGAL_FAILED}: {e}")

    try:
        out_path = os.path.join(shim.ready_genome_files_dir, f"{gd.id}.faa")
        num = _filter_and_rename_fasta(out_path, raw_aa_path, gd.id, 99999)
        if num == 0:
            remove_file_if_exists(out_path)
            raise TargetGenomeError(MISSED_NO_PROTEINS)
        return out_path
    finally:
        remove_file_if_exists(raw_aa_path)
