from gtotree.utils.messaging import report_processing_stage
# from gtotree.utils.seqs

def update_headers(run_data):

    if not run_data.updating_headers:
        return run_data

    report_processing_stage("updating-headers")

    ### functions to update run_data.mapping_dict
        # keep a master initial list of all entries that were renamed by user-provided mapping file
            # so we know not to change those with taxonomic info
    # run_data = update_with_ncbi_tax(run_data)
        # only if no mapping_dict info from a specific mapping file for a given entry
    # run_data = update_with_gtdb_tax(run_data)
        # only if no mapping_dict info from a specific mapping file for a given entry

    ### at this point the run_data.mapping_dict should have all entries that need to be changed
    # run_data = create_modified_alignment_file(run_data)
        # write out to, e.g., os.path.join(run_data.output_dir, "aligned-SCGs-mod-names.faa")

    return run_data
