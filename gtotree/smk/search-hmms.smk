import os
from gtotree.utils.general import (read_run_data,
                                   write_run_data)
from gtotree.utils.hmms.hmm_searching import (run_hmm_search,
                                             parse_hmmer_results)

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

genome_dict = {gd.id: gd for gd in run_data.get_all_input_genome_for_hmm_search()}
genome_ids = list(genome_dict.keys())

rule all:
    input:
        expand(f"{run_data.hmm_results_dir}/{{ID}}.done", ID=genome_ids)
    run:
        for ID in genome_ids:
            genome = genome_dict[ID]
            status_path = f"{run_data.hmm_results_dir}/{ID}.done"

            with open(status_path, 'r') as f:
                for line in f:
                    ID, hmm_search_failed = line.strip().split('\t')

                    if int(hmm_search_failed):
                        genome.mark_hmm_search_failed()
                        genome.mark_removed()
                    else:
                        genome.mark_hmm_search_done()

        write_run_data(run_data)


rule search_hmms:
    output:
        f"{run_data.hmm_results_dir}/{{ID}}.done"
    run:
        genome = genome_dict[wildcards.ID]
        AA_path = genome.final_AA_path

        hmm_search_failed = run_hmm_search(wildcards.ID, run_data, AA_path)
        dict_of_hit_counts, dict_of_hit_gene_ids = parse_hmmer_results(f"{run_data.hmm_results_dir}/{wildcards.ID}-hmm-hits.txt",
                                                                       run_data)

        ## esl-sfetch index
        ## then iterate over dict_of_hit_gene_ids,
            ## if gene_id present
                ## esl-sfetch again to get the seq
                ## rename header to genome_id
                ## write to file (i think it's own for the genome right now, then have to merge in "all")
                    ## maybe put in subdirectory for each genome to keep the number of files per dir down a bit
        print(dict_of_hit_gene_ids)

        with open(output[0], 'w') as f:
            f.write(f"{wildcards.ID}\t{int(hmm_search_failed)}\n")

