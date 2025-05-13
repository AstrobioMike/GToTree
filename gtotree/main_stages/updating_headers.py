from gtotree.utils.messaging import report_processing_stage
from gtotree.utils.gtdb.handle_gtdb_tax_info import update_mapping_dict_with_gtdb_tax_info
from gtotree.utils.ncbi.handle_ncbi_tax_info import update_mapping_dict_with_ncbi_tax_info
from gtotree.utils.seqs import swap_labels_in_alignment
from gtotree.utils.general import write_run_data

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

        run_data = swap_labels_in_alignment(run_data)
        run_data.headers_updated = True
        write_run_data(run_data)

    return run_data
