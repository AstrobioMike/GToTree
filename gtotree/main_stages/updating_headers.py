from gtotree.utils.misc.messaging import report_processing_stage, report_early_exit
from gtotree.utils.gtdb.handle_gtdb_tax_info import update_mapping_dict_with_gtdb_tax_info
from gtotree.utils.ncbi.handle_ncbi_tax_info import update_mapping_dict_with_ncbi_tax_info
from gtotree.utils.misc.seqs import swap_labels_in_alignment
from gtotree.utils.misc.general import write_run_data

def update_headers(args, run_data):

    if not run_data.updating_headers:
        run_data.final_alignment_path = run_data.concatenated_alignment_path
        return run_data

    report_processing_stage("updating-headers", run_data)

    if not run_data.headers_updated:

        if args.add_gtdb_tax:
            run_data = update_mapping_dict_with_gtdb_tax_info(args, run_data)

        if args.add_ncbi_tax:
            run_data = update_mapping_dict_with_ncbi_tax_info(args, run_data)

        try:
            run_data = swap_labels_in_alignment(run_data)
        except FileNotFoundError as e:
            report_early_exit(
                run_data,
                "We were about to swap in the more informative headers for some "
                f"reason. If the output directory was modified since the previous "
                "run, your best bet is to start fresh with `-F`.")

        run_data.headers_updated = True
        write_run_data(run_data)

    return run_data
