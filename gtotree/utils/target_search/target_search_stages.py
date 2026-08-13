"""
The work phases behind `gtt search-pfams` / `gtt search-kos`
"""

import os
import tempfile
from tqdm import tqdm  # type: ignore
from gtotree.utils.misc import phase_stats
from gtotree.utils.misc.general import (run_pooled_stage,
                                   write_run_data,
                                   resolve_input_genomes as shared_resolve_input_genomes,
                                   SOURCE_ACCESSION, SOURCE_GENBANK, SOURCE_FASTA,
                                   SOURCE_AMINO_ACID,
                                   GTT_PROGRESS_BAR_FORMAT_INDENTED,
                                   GTT_PROGRESS_BAR_FORMAT_INDENTED_6,
                                   GTT_PROGRESS_SMOOTHING)
from gtotree.utils.misc.messaging import (report_message, color_text, spinner,
                                          REMOVED_GENOMES_FILENAME)
from gtotree.utils.misc.summary_info import write_removed_genomes_report
from gtotree.utils.misc.stages import (GenomeRemovalStage,
                                       GENOME_REMOVAL_STAGE_ORDER)
from gtotree.utils.hmms.hmm_searching_engine import press_profiles
from gtotree.main_stages.processing_genomes import (SearchPlan,
                                                    _fused,
                                                    genomes_needing_processing)
from gtotree.utils.misc.processing_genomes import (
    build_base_link_map,
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


def _removals_pointer(run_data):
    return f"{run_data.run_files_dir_rel}/{REMOVED_GENOMES_FILENAME}"


SEARCH_PHASE_REMOVAL_STAGES = tuple(
    stage for stage in GENOME_REMOVAL_STAGE_ORDER
    if stage != GenomeRemovalStage.NCBI_LOOKUP)


################################################################################
# phase 1: resolving the input genome set
################################################################################

def resolve_input_genomes(args, run_data):
    """
    Fold in any `-w` selections, then report the input genome set

    The body is `general.resolve_input_genomes`, shared with `gtt gen-scg-hmms` so all
    three programs render this phase identically. Only the exception type differs.

    Returns (run_data, selections), where `selections` is the list of
    RefGenomeSelection objects `-w` produced (empty when `-w` wasn't used).
    """
    return shared_resolve_input_genomes(args, run_data, TargetSearchError)


def lookup_ncbi_accessions(run_data):
    """
    Resolve accessions against the NCBI assembly summary

    This supplies the download links the processing stage needs, and marks any
    accession NCBI no longer lists as removed, so it lands in the summary table with a
    reason rather than silently failing to download later.
    """
    from gtotree.utils.ncbi.parse_ncbi_assembly_summary import parse_assembly_summary
    from gtotree.utils.ncbi.get_ncbi_assembly_data import get_ncbi_assembly_summary_tab

    if not run_data.ncbi_accs:
        return run_data

    with spinner("Looking up accessions in the NCBI assembly data...",
                 "Looked up accessions at NCBI"):
        run_data = parse_assembly_summary(get_ncbi_assembly_summary_tab(), run_data)

    not_found = run_data.get_ncbi_accs_not_found()
    if not_found:
        report_message(
            f"{len(not_found):,} accession(s) not found at NCBI. Reported in:",
            "yellow", ii="      ", si="      ")
        print(f"        {color_text(_removals_pointer(run_data), 'yellow')}")

        remaining = len([gd for gd in run_data.all_input_genomes if not gd.removed])
        if not remaining:
            raise TargetSearchError(
                "None of the input accessions were found at NCBI, so there's nothing "
                "left to search.")

    return run_data


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
    Resolve the requested target IDs to searchable profiles

    For Pfams this is a streaming pass over the master Pfam-A.hmm pulling out the
    wanted profiles; for KOs it's a lookup in ko_list plus copying the matching HMMs.
    Both are the same functions the main driver calls, and both write their own
    "requested vs found" bookkeeping into the output directory.
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

    print(f"      {len(found):,} of {total:,} requested "
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
            f"`{spec.targets_flag}`, they should look like "
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
    return SearchPlan(
        do_pfam=(spec.plan_flag == "do_pfam"),
        do_ko=(spec.plan_flag == "do_ko"),
        keep_genome_files=bool(getattr(args, "debug", False)),
        do_scg=False,
    )


# GenomeData.source -> (preprocessing worker, status applier). The accession entry is
# added per run, since its worker closes over that run's base-link map
_SOURCE_WORKERS = {
    SOURCE_GENBANK: (_process_one_genbank_file, _apply_genbank_status),
    SOURCE_FASTA: (_process_one_fasta_file, _apply_fasta_status),
    SOURCE_AMINO_ACID: (_process_one_amino_acid_file, _apply_amino_acid_status),
}


def _dispatching_worker_pair(run_data):
    """
    One preprocessing worker/applier pair covering all four input sources

    `run_pooled_stage` takes a single worker, so the per-source dispatch happens inside
    it rather than by running four pools. The pair goes through `_fused` exactly like
    the main driver's per-source pairs do, so each genome is still searched and dropped
    in the same worker that produced it.

    The base-link map is built once here rather than per genome.
    """
    workers = dict(_SOURCE_WORKERS)

    if run_data.ncbi_accs:
        base_link_map = build_base_link_map(run_data)

        def _ncbi(gd, rd):
            return _process_one_ncbi_accession(gd, rd, base_link_map)

        workers[SOURCE_ACCESSION] = (_ncbi, _apply_ncbi_accession_status)

    def preprocess(gd, rd):
        return workers[gd.source][0](gd, rd)

    def apply_status(gd, status, rd):
        workers[gd.source][1](gd, status, rd)

    return preprocess, apply_status


def phase_search_genomes(args, run_data, spec, plan):
    """
    Preprocess and search every genome that isn't already done, in one pool with one
    progress bar

    Genomes already carrying both `processing_done` and this target type's search flag
    are skipped, which is what makes `-R` resume at genome granularity without any
    extra bookkeeping. Accessions that weren't found at NCBI were marked removed back
    in phase 1, so they're excluded here by `genomes_needing_processing` rather than
    needing a filter of their own.
    """
    phase_stats.begin("processing and searching genomes")

    press_dir_cm = (tempfile.TemporaryDirectory(prefix="gtt-press-")
                    if spec.presses_profiles else _NullContext())

    with press_dir_cm as press_dir:
        if spec.presses_profiles:
            hmm_path = getattr(run_data, spec.combined_hmm_attr)
            with spinner(f"Preparing {spec.target_label} profiles for searching...",
                         f"Prepared {spec.target_label} profiles"):
                plan.pressed_pfam_base = press_profiles(
                    hmm_path, press_dir, "target-profiles")

        to_process = genomes_needing_processing(run_data.all_input_genomes, plan)

        # only genomes still in the run count towards "already done"; ones dropped at
        # the NCBI lookup are a different thing entirely and were reported in phase 1
        alive = [gd for gd in run_data.all_input_genomes if not gd.removed]
        already_done = len(alive) - len(to_process)

        num_targets = len(spec.found_targets(run_data))

        if not to_process:
            print(f"\n      All {len(alive):,} genome(s) were already processed and "
                  "searched in a previous run")
            return run_data

        print(f"\n      Processing and searching {len(to_process):,} genome(s) for "
              f"{num_targets:,} {spec.target_label} target(s):")
        if already_done:
            print(f"        ({already_done:,} already done in a previous run)")

        preprocess, apply_status = _dispatching_worker_pair(run_data)
        worker, apply_result = _fused(preprocess, apply_status, plan)

        run_data = run_pooled_stage(to_process, worker, apply_result, args, run_data,
                                    bar_format=GTT_PROGRESS_BAR_FORMAT_INDENTED_6)

    mark_failed_searches_removed(run_data, spec)

    write_removed_genomes_report(run_data)
    write_run_data(run_data)

    searched = count_searched_genomes(run_data, spec)
    failed_here = run_data.genomes_removed_at(*SEARCH_PHASE_REMOVAL_STAGES)

    # print(f"\n      {color_text(f'{searched:,} genome(s) searched', 'green')}")

    if failed_here:
        report_message(
            f"{len(failed_here):,} genome(s) failed the search phase. Reported in:",
            "yellow", ii="      ", si="      ")
        print(f"        {color_text(_removals_pointer(run_data), 'yellow')}")

    if not searched:
        raise TargetSearchError(
            "None of the input genomes made it through to a completed search. The "
            f"reason for each is in {_removals_pointer(run_data)}.")

    return run_data


class _NullContext:
    """Stand-in for the press dir when this target type doesn't press profiles."""

    def __enter__(self):
        return None

    def __exit__(self, *exc):
        return False


################################################################################
# phase 4: outputs
################################################################################

def phase_write_outputs(run_data, spec, out_dir):
    """
    Build the combined outputs from the per-genome artifacts

    Both helpers are the main driver's, unmodified: they read the per-genome result
    files off disk rather than any in-memory accumulation, which is why they're
    idempotent and complete regardless of how much of the work happened in this
    invocation versus an earlier interrupted one.

    Returns (summary_path, num_targets_with_hits).
    """
    from gtotree.utils.target_search import target_search_outputs as outputs

    run_data.update_all_input_genomes()

    # with spinner("Writing the hit-counts table...", "Wrote the hit-counts table"):
    #     spec.write_counts_table(run_data)

    spec.write_counts_table(run_data)

    found = spec.found_targets(run_data)
    hit_seqs_dir = os.path.join(spec.results_dir(run_data), spec.hit_seqs_subdir)

    # print("      Combining per-genome hit sequences:")
    with tqdm(total=len(found), bar_format=GTT_PROGRESS_BAR_FORMAT_INDENTED_6,
              ncols=76, smoothing=GTT_PROGRESS_SMOOTHING) as pbar:
        # combine_hits takes the whole target list, so it's driven one target at a
        # time here purely to give the bar something to advance on
        for target in found:
            spec.combine_hits([target], spec.tmp_results_dir(run_data), hit_seqs_dir)
            pbar.update(1)

    print()

    targets_with_hits = _count_and_prune_hit_seqs(hit_seqs_dir)

    # with spinner("Writing the genome summary table...", "Wrote the genome summary table"):
    #     summary_path = outputs.write_genomes_summary(out_dir, run_data, spec)

    summary_path = outputs.write_genomes_summary(out_dir, run_data, spec)

    write_run_data(run_data)

    return summary_path, targets_with_hits


def _count_and_prune_hit_seqs(hit_seqs_dir):
    """
    How many targets got a combined hit-seqs file, removing the directory when that's
    none of them

    `combine_hits` writes nothing for a target with no hits, so an empty directory here
    means not one genome hit not one target. Leaving it behind reads like the run
    produced sequences that then went missing, and the finish banner would point at it.
    """
    if not os.path.isdir(hit_seqs_dir):
        return 0

    contents = os.listdir(hit_seqs_dir)

    if not contents:
        try:
            os.rmdir(hit_seqs_dir)
        except OSError:
            pass
        return 0

    return len(contents)


def mark_failed_searches_removed(run_data, spec):
    """
    Drop genomes whose search failed at the search stage
    """
    for gd in run_data.all_input_genomes:
        if gd.removed:
            continue
        if getattr(gd, spec.search_failed_flag, False):
            gd.mark_removed(f"{spec.target_label} search failed",
                            GenomeRemovalStage.HMM_SEARCH)

    return run_data


def count_searched_genomes(run_data, spec):
    """
    How many genomes actually made it all the way through to a completed search
    """
    return len([gd for gd in run_data.all_input_genomes
                if getattr(gd, spec.search_done_flag, False)
                and not getattr(gd, spec.search_failed_flag, False)
                and not gd.removed])
