from gtotree.utils.messaging import report_processing_stage
from gtotree.utils.seqs import (concatenate_alignments,
                                gen_partitions_file,
                                get_alignment_length)

def concatenate_SCG_sets(run_data):

    report_processing_stage("concatenate-SCG-sets", run_data)

    run_data, dict_of_genomes, SCG_IDs = concatenate_alignments(run_data)

    run_data.final_alignment_length = get_alignment_length(run_data.concatenated_alignment_path)

    gen_partitions_file(run_data, SCG_IDs, dict_of_genomes)

    return run_data
