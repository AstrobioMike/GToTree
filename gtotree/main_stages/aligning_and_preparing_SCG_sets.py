import os
from gtotree.utils.general import (write_run_data,
                                   run_pooled_stage,
                                   remove_file_if_exists,
                                   check_file_exists_and_not_empty)
from gtotree.utils.messaging import (report_processing_stage,
                                     report_no_SCGs_remaining)
from gtotree.utils.seqs import (copy_gene_alignments,
                                run_muscle,
                                run_trimal,
                                add_needed_gap_seqs,
                                fasta_has_single_record)


def align_and_prepare_SCG_sets(args, run_data):

    report_processing_stage("align-and-prepare-SCG-sets", run_data)

    if not run_data.all_SCG_sets_aligned:

        scgs_to_align = run_data.get_all_SCG_targets_remaining()

        if len(scgs_to_align) > 0:
            run_data = run_pooled_stage(scgs_to_align,
                                        _align_one_SCG,
                                        _apply_alignment_status,
                                        args, run_data)
            run_data.all_SCG_sets_aligned = True
            write_run_data(run_data)

        else:
            report_no_SCGs_remaining(run_data)

        if args.keep_gene_alignments:
            copy_gene_alignments(run_data)

    return run_data


def _scg_paths(scg_id, run_data):
    """All per-SCG paths for the align/trim/finalize steps, in one place."""
    d = run_data.found_SCG_seqs_dir
    ext = run_data.general_ext
    return {
        "inpath":     f"{d}/{scg_id}-genome-filtered{ext}",
        "aligned":    f"{d}/{scg_id}-aligned{ext}",
        "align_log":  f"{d}/{scg_id}-align.log",
        "trimmed":    f"{d}/{scg_id}-trimmed{ext}",
        "trimal_log": f"{d}/{scg_id}-trimmal.log",
        "final":      f"{d}/{scg_id}-final{ext}",
    }


def _align_one_SCG(scg, run_data):
    """
    Worker: align + trim + gap-pad one SCG set. Touches only this SCG's
    files and returns a status dict; all SCGset mutation happens on the main thread
    in _apply_alignment_status.

    Restart-safety:
      - a valid (non-empty) '-final' from a prior run means this SCG is done ->
        return skipped, so a re-run resumes rather than realigning everything
      - before doing any work, stale '.part' temp files and stale intermediate
        products ('-aligned', '-trimmed') from a prior *killed* run are removed so
        muscle/trimal never consume a truncated intermediate
      - '-final' is written atomically (build '-final.part', os.replace on success;
        remove the '.part' on any failure) so a killed run never leaves a truncated
        '-final' that a later skip-check would wrongly trust
    """
    paths = _scg_paths(scg.id, run_data)
    final_path = paths["final"]
    final_part = final_path + ".part"

    # a valid final output already exists -> this SCG was completed on a prior run
    if check_file_exists_and_not_empty(final_path):
        return {"skipped": True, "align_failed": False, "trimal_failed": False}

    # clean up anything left behind by a prior interrupted run before starting:
    # the leftover atomic temp, plus intermediates we're about to regenerate
    # (check_file_exists_and_not_empty already drops a zero-byte -final for us)
    remove_file_if_exists(final_part)
    remove_file_if_exists(paths["aligned"])
    remove_file_if_exists(paths["trimmed"])

    try:
        if fasta_has_single_record(paths["inpath"]):
            align_failed = False
            trimal_failed = False
            trimmed_source = paths["inpath"]
            with open(paths["align_log"], "w") as f:
                f.write("Single sequence in SCG set; skipping alignment.\n")
        else:
            align_failed = run_muscle(scg.id, run_data, paths["inpath"],
                                      paths["aligned"], paths["align_log"])
            trimal_failed = run_trimal(paths["aligned"], paths["trimmed"],
                                       paths["trimal_log"])
            trimmed_source = paths["trimmed"]

        # only finalize when both upstream steps succeeded; otherwise leave no
        # -final so the failure is visible and a re-run retries this SCG
        if not align_failed and not trimal_failed:
            add_needed_gap_seqs(run_data, trimmed_source, final_part)
            os.replace(final_part, final_path)
        else:
            remove_file_if_exists(final_part)

        return {
            "skipped": False,
            "align_failed": bool(align_failed),
            "trimal_failed": bool(trimal_failed),
        }

    except BaseException:
        # never leave a half-written final behind, and never take down the pool;
        # surface as an alignment failure so bookkeeping marks the SCG removed
        remove_file_if_exists(final_part)
        return {"skipped": False, "align_failed": True, "trimal_failed": False}


def _apply_alignment_status(scg, status, run_data):
    """
    Apply one worker's result to its SCGset. Main thread only. A skipped SCG (its
    '-final' already existed) is treated as fully aligned/trimmed/ready.
    """
    if status.get("skipped"):
        scg.aligned = True
        scg.trimmed = True
        scg.ready_for_cat = True
        return

    if status.get("align_failed"):
        scg.mark_removed("alignment failed")
        scg.aligned = False
    elif status.get("trimal_failed"):
        scg.mark_removed("trimal failed")
        scg.trimmed = False
    else:
        scg.aligned = True
        scg.trimmed = True
        scg.ready_for_cat = True
