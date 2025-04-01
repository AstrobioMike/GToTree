from gtotree.utils.messaging import report_processing_stage
from gtotree.utils.seqs import concatenate_alignments, gen_partitions_file

def concatenate_SCG_sets(run_data):

    report_processing_stage("concatenate-SCG-sets")

    dict_of_genomes, SCG_IDs = concatenate_alignments(run_data)

    gen_partitions_file(run_data, SCG_IDs, dict_of_genomes)

    return run_data
