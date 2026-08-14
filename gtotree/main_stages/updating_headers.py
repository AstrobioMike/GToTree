from gtotree.utils.misc.messaging import (report_processing_stage, report_early_exit,
                                          report_message)
from gtotree.utils.gtdb.handle_gtdb_tax_info import update_mapping_dict_with_gtdb_tax_info
from gtotree.utils.ncbi.handle_ncbi_tax_info import update_mapping_dict_with_ncbi_tax_info
from gtotree.utils.misc.seqs import (swap_labels_in_alignment,
                                     usable_concatenated_alignment_path)
from gtotree.utils.misc.general import write_run_data
from gtotree.utils.misc.stages import PipelineStage

def update_headers(args, run_data):

    if not run_data.updating_headers:
        run_data.final_alignment_path = run_data.concatenated_alignment_path
        return run_data

    report_processing_stage("updating-headers", run_data)

    if not run_data.stage_is_complete(PipelineStage.UPDATE_HEADERS):

        run_data.header_update_error = ""

        try:
            if args.add_gtdb_tax:
                run_data = update_mapping_dict_with_gtdb_tax_info(args, run_data)

            if args.add_ncbi_tax:
                run_data = update_mapping_dict_with_ncbi_tax_info(args, run_data)

            run_data = swap_labels_in_alignment(run_data)

        except Exception as err:
            return _carry_on_with_original_labels(run_data, err)

        run_data.mark_stage_complete(PipelineStage.UPDATE_HEADERS)
        write_run_data(run_data)

    return run_data


def _carry_on_with_original_labels(run_data, err):
    """
    Salvage the run if possible if header swapping fails

    Swapping in fuller labels is a convenience layered on top of a finished alignment,
    and that failing doesn't necessarily mean we want to throw away the whole run
    thar already aligned and concatenated things.

    The stage is deliberately NOT marked complete, so a later `-R` retries the swap
    rather than inheriting the fallback.
    """
    fallback_path = usable_concatenated_alignment_path(run_data)

    if fallback_path is None:
        report_early_exit(
            run_data,
            "We couldn't swap in the more informative headers because "
            f"{err}. If the output directory was modified since the previous "
            "run, your best bet is to start fresh with `-F`.")

    run_data.header_update_error = str(err)
    run_data.final_alignment_path = fallback_path

    report_message(
        "We couldn't swap in the more informative headers because "
        f"{err}. Rather than lose the run over it, we're carrying on with the "
        "original genome IDs, so the final alignment and tree will hold those "
        "instead of the swapped labels.", "yellow")

    write_run_data(run_data)

    return run_data
