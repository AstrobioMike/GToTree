from gtotree.utils.messaging import report_processing_stage
from gtotree.utils.gtdb.handle_gtdb_tax_info import update_mapping_dict_with_gtdb_tax_info
from gtotree.utils.ncbi.handle_ncbi_tax_info import update_mapping_dict_with_ncbi_tax_info
# from gtotree.utils.seqs

def update_headers(args, run_data):

    if not run_data.updating_headers:
        return run_data

    report_processing_stage("updating-headers")

    if args.add_gtdb_tax:
        run_data = update_mapping_dict_with_gtdb_tax_info(args, run_data)

    if args.add_ncbi_tax:
        run_data = update_mapping_dict_with_ncbi_tax_info(args, run_data)


    ### at this point the run_data.mapping_dict should have all entries that need to be changed
    # run_data = create_modified_alignment_file(run_data)
        # write out to, e.g., os.path.join(run_data.output_dir, "aligned-SCGs-mod-names.faa")

    return run_data
