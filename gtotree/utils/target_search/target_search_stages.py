"""
The work phases behind `gtt search-pfams` / `gtt search-kos`.

The per-genome heavy lifting is not reimplemented here. `main_stages.processing_genomes`
already owns a fused worker that, for one genome, preprocesses it, immediately searches
the resulting FASTA while it's still on disk, and then drops the FASTA -- which is what
keeps peak disk usage flat across a large run. That worker dispatches on
`SearchPlan.do_pfam` / `do_ko`, so these subcommands get it by constructing a plan with
`do_scg=False` and the relevant target flag set, rather than by growing a second copy.

What *is* here is the orchestration and reporting layer: the phase structure, spinners,
and progress bars follow `gtt gen-scg-hmms` rather than the main GToTree driver, whose
`report_*_update` blocks assume an SCG set, a tree, and a four-genome floor.
"""

import os
import tempfile

from tqdm import tqdm  # type: ignore

from gtotree.utils.misc import phase_stats
from gtotree.utils.misc.general import (run_pooled_stage,
                                   write_run_data,
                                   GTT_PROGRESS_BAR_FORMAT_INDENTED,
                                   GTT_PROGRESS_SMOOTHING)
from gtotree.utils.misc.messaging import report_message, color_text, spinner
from gtotree.utils.hmms.hmm_searching_engine import press_profiles
from gtotree.main_stages.processing_genomes import (SearchPlan,
                                                    _fused,
                                                    genomes_needing_processing)
from gtotree.utils.misc.processing_genomes import (
    build_base_link_map,
    capture_ncbi_failed_downloads,
    capture_failed_genbank_files,
    capture_failed_fasta_files,
    capture_failed_amino_acid_files,
    _process_one_ncbi_accession,
    _apply_ncbi_accession_status,
    _process_one_genbank_file,
    _apply_genbank_status,
    _process_one_fasta_file,
    _apply_fasta_status,
    _process_one_amino_acid_file,
    _apply_amino_acid_status,
)
from gtotree.utils.target_search.target_search_setup import TargetSearchError


################################################################################
# phase 1: input genomes
################################################################################

def phase_resolve_genomes(args, run_data):
    """
    Report the input genomes already parsed into run_data, and fold in any `-w`
    selection.

    Returns (run_data, selection) where `selection` is the RefGenomeSelection when `-w`
    was used, else None. The caller needs it for the run-state fingerprint.
    """
    from gtotree.utils.taxonomy.wanted_ref_tax import (resolve_wanted_ref_tax_accessions,
                                                       describe_source_version)

    counts = [
        (args.ncbi_accessions, len(run_data.get_user_provided_ncbi_accs()),
         "NCBI accession(s)"),
        (args.genbank_files, len(run_data.genbank_files), "GenBank file(s)"),
        (args.fasta_files, len(run_data.fasta_files), "nucleotide fasta file(s)"),
        (args.amino_acid_files, len(run_data.amino_acid_files),
         "amino-acid fasta file(s)"),
    ]

    for provided, count, label in counts:
        if provided:
            print(f"        {count:,} {label} read from {provided}")

    selection = None

    if args.wanted_ref_tax:
        source_desc = describe_source_version(args.source)
        if source_desc:
            print(f"\n      Genome source being used: "
                  f"{color_text(source_desc, 'green')}\n")

        with spinner(f"Selecting reference genomes for '{args.wanted_ref_tax}'...",
                     "Selected reference genomes"):
            accessions, selection = resolve_wanted_ref_tax_accessions(
                args.source, args.wanted_ref_tax,
                target_rank=args.target_rank,
                derep_rank=args.derep_rank)

        added = run_data.merge_wanted_ref_tax_accessions(accessions)

        detail = f"        {len(accessions):,} genome(s) selected"
        if selection.resolved_rank:
            detail += f" ({selection.canonical} at rank {selection.resolved_rank})"
        print(detail)
        if selection.effective_derep_rank:
            print("        dereplicated to one genome per "
                  f"{selection.effective_derep_rank}")
        if added != len(accessions):
            print(f"        {len(accessions) - added:,} already present from `-a`")

        for warning in selection.warnings:
            report_message(warning, "orange", ii="        ", si="        ")

    run_data.update_all_input_genomes()
    total = len(run_data.all_input_genomes)

    if not total:
        raise TargetSearchError("No input genomes were resolved to work with.")

    print(f"\n      {color_text(f'{total:,} total input genome(s)', 'green')}")

    return run_data, selection


################################################################################
# phase 2: target collection
################################################################################

def target_stage_is_reusable(run_data, spec, out_dir):
    """
    True when a resumed run can trust the previously collected target set.

    The done-flag alone isn't enough: it lives in run-data.json, but the profiles it
    describes live on disk in the output directory, and the two can drift (an
    interrupted write, a user tidying up between runs). So the flag has to be backed by
    the artifacts actually being there.
    """
    if not getattr(run_data, spec.searching_done_attr, False):
        return False
    if not spec.found_targets(run_data):
        return False
    for rel in spec.target_stage_artifacts:
        if not os.path.exists(os.path.join(out_dir, rel)):
            return False
    return True


def phase_collect_targets(run_data, spec, out_dir, resuming=False):
    """
    Resolve the requested target IDs to searchable profiles.

    For Pfams this is a streaming pass over the master Pfam-A.hmm pulling out the
    wanted profiles; for KOs it's a lookup in ko_list plus copying the matching HMMs.
    Both are the same functions the main driver calls, and both write their own
    "requested vs found" bookkeeping into the output directory.

    Reusing this on resume matters: the Pfam pass in particular is a multi-minute scan
    over a ~2 GB file.
    """
    if resuming and target_stage_is_reusable(run_data, spec, out_dir):
        found = spec.found_targets(run_data)
        with spinner(f"Reusing previously collected {spec.target_label} targets...",
                     f"Reused {len(found):,} {spec.target_label} target(s)"):
            pass
    else:
        setattr(run_data, spec.searching_done_attr, False)
        run_data = spec.collect_targets(run_data)
        setattr(run_data, spec.searching_done_attr, True)

    found = spec.found_targets(run_data)
    failed = spec.failed_targets(run_data)
    total = spec.total_targets(run_data)

    print(f"        {len(found):,} of {total:,} requested "
          f"{spec.target_label_plural} found and ready to search")

    if failed:
        report_message(
            f"{len(failed):,} requested {spec.target_label} target(s) couldn't be "
            f"found, and are reported in:", "yellow", ii="        ", si="        ")
        print(f"          {color_text(spec.failed_targets_filename, 'yellow')}")

    if not found:
        raise TargetSearchError(
            f"None of the requested {spec.target_label_plural} could be found, so "
            "there's nothing to search the genomes for. Check the IDs in "
            f"`{spec.targets_flag}` -- they should look like "
            f"'{spec.example_target}'.")

    return run_data


################################################################################
# phase 3: processing + searching
################################################################################

def build_plan(args, spec):
    """
    A SearchPlan that turns on exactly this subcommand's search and nothing else.

    `do_scg=False` is what makes the shared fused worker usable here: there's no SCG
    HMM set to press or search, and no tree downstream that would consume one.
    """
    plan = SearchPlan(
        do_pfam=(spec.plan_flag == "do_pfam"),
        do_ko=(spec.plan_flag == "do_ko"),
        keep_genome_files=bool(getattr(args, "debug", False)),
        do_scg=False,
    )
    return plan


def phase_search_genomes(args, run_data, spec, plan):
    """
    Preprocess and search every genome that isn't already done.

    Each source gets its own labelled progress bar. Genomes already carrying both
    `processing_done` and this target type's search flag are skipped, which is what
    makes `-R` resume at genome granularity without any extra bookkeeping.
    """
    press_dir_cm = (tempfile.TemporaryDirectory(prefix="gtt-press-")
                    if spec.presses_profiles else _NullContext())

    with press_dir_cm as press_dir:
        if spec.presses_profiles:
            hmm_path = getattr(run_data, spec.combined_hmm_attr)
            with spinner(f"Preparing {spec.target_label} profiles for searching...",
                         f"Prepared {spec.target_label} profiles"):
                plan.pressed_pfam_base = press_profiles(
                    hmm_path, press_dir, "target-profiles")

        run_data = _run_source(args, run_data, plan, spec,
                               source_list=run_data.ncbi_accs,
                               label="NCBI accessions",
                               phase="ncbi accessions",
                               prepare=_prepare_ncbi,
                               capture=capture_ncbi_failed_downloads)

        run_data = _run_source(args, run_data, plan, spec,
                               source_list=run_data.genbank_files,
                               label="GenBank files",
                               phase="genbank files",
                               worker_pair=(_process_one_genbank_file,
                                            _apply_genbank_status),
                               capture=capture_failed_genbank_files)

        run_data = _run_source(args, run_data, plan, spec,
                               source_list=run_data.fasta_files,
                               label="fasta files",
                               phase="fasta files",
                               worker_pair=(_process_one_fasta_file,
                                            _apply_fasta_status),
                               capture=capture_failed_fasta_files)

        run_data = _run_source(args, run_data, plan, spec,
                               source_list=run_data.amino_acid_files,
                               label="amino-acid files",
                               phase="amino-acid files",
                               worker_pair=(_process_one_amino_acid_file,
                                            _apply_amino_acid_status),
                               capture=capture_failed_amino_acid_files)

    return run_data


class _NullContext:
    """Stand-in for the press dir when this target type doesn't press profiles."""

    def __enter__(self):
        return None

    def __exit__(self, *exc):
        return False


def _prepare_ncbi(args, run_data):
    """
    Resolve accessions against the NCBI assembly summary before downloading.

    This both supplies the download links and marks any accession NCBI no longer lists
    as removed, so it shows up in the summary table with a reason rather than silently
    failing to download.
    """
    from gtotree.utils.ncbi.parse_ncbi_assembly_summary import parse_assembly_summary
    from gtotree.utils.ncbi.get_ncbi_assembly_data import get_ncbi_assembly_summary_tab

    with spinner("Looking up accessions in the NCBI assembly data...",
                 "Looked up accessions at NCBI"):
        run_data = parse_assembly_summary(get_ncbi_assembly_summary_tab(), run_data)

    base_link_map = build_base_link_map(run_data)

    def preprocess(acc_gd, rd):
        return _process_one_ncbi_accession(acc_gd, rd, base_link_map)

    return run_data, (preprocess, _apply_ncbi_accession_status), _only_found


def _only_found(genomes):
    return [gd for gd in genomes if gd.acc_was_found]


def _run_source(args, run_data, plan, spec, source_list, label, phase,
                worker_pair=None, prepare=None, capture=None, ):
    """
    Run one input source through the fused preprocess-and-search stage.
    """
    if not source_list:
        return run_data

    phase_stats.begin(f"searching: {phase}")

    extra_filter = None
    if prepare is not None:
        run_data, worker_pair, extra_filter = prepare(args, run_data)

    to_process = genomes_needing_processing(source_list, plan)
    if extra_filter is not None:
        to_process = extra_filter(to_process)

    already_done = len(source_list) - len(to_process)

    if not to_process:
        print(f"\n      All {len(source_list):,} {label} were already done")
        return run_data

    print(f"\n      Processing and searching {len(to_process):,} {label}:")
    if already_done:
        print(f"        ({already_done:,} already done in a previous run)")

    preprocess, apply_status = worker_pair
    worker, apply_result = _fused(preprocess, apply_status, plan)

    run_data = run_pooled_stage(to_process, worker, apply_result, args, run_data,
                                bar_format=GTT_PROGRESS_BAR_FORMAT_INDENTED)

    if capture is not None:
        capture(run_data)

    write_run_data(run_data)

    return run_data


################################################################################
# phase 4: outputs
################################################################################

def phase_write_outputs(run_data, spec, out_dir):
    """
    Build the combined outputs from the per-genome artifacts.

    Both helpers are the main driver's, unmodified: they read the per-genome result
    files off disk rather than any in-memory accumulation, which is why they're
    idempotent and complete regardless of how much of the work happened in this
    invocation versus an earlier interrupted one.
    """
    from gtotree.utils.target_search import target_search_outputs as outputs

    run_data.update_all_input_genomes()

    with spinner("Writing the hit-counts table...", "Wrote the hit-counts table"):
        spec.write_counts_table(run_data)

    found = spec.found_targets(run_data)
    hit_seqs_dir = os.path.join(spec.results_dir(run_data), spec.hit_seqs_subdir)

    print("\n      Combining per-genome hit sequences:")
    with tqdm(total=len(found), bar_format=GTT_PROGRESS_BAR_FORMAT_INDENTED,
              ncols=76, smoothing=GTT_PROGRESS_SMOOTHING) as pbar:
        # combine_hits takes the whole target list, so it's driven one target at a
        # time here purely to give the bar something to advance on
        for target in found:
            spec.combine_hits([target], spec.tmp_results_dir(run_data), hit_seqs_dir)
            pbar.update(1)

    with spinner("Writing the genome summary table...", "Wrote the genome summary table"):
        summary_path = outputs.write_genomes_summary(out_dir, run_data, spec)

    write_run_data(run_data)

    return summary_path


def count_searched_genomes(run_data, spec):
    """
    How many genomes actually made it all the way through to a completed search.
    """
    return len([gd for gd in run_data.all_input_genomes
                if getattr(gd, spec.search_done_flag, False) and not gd.removed])
