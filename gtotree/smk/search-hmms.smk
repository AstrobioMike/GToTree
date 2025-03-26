import os
from gtotree.utils.general import (read_run_data,
                                   write_run_data,
                                   run_hmm_search)

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

genome_dict = {gd.id: gd for gd in run_data.get_all_remaining_input_genomes()}
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
                    ID, status = line.strip().split('\t')

                    if not int(status):
                        genome.mark_hmm_search_failed()
                        genome.mark_removed()

rule search_hmms:
    output:
        f"{run_data.hmm_results_dir}/{{ID}}.done"
    run:
        genome = genome_dict[wildcards.ID]
        AA_path = genome.final_AA_path

        done = run_hmm_search(wildcards.ID, run_data, AA_path)

        with open(output[0], 'w') as f:
            f.write(f"{wildcards.ID}\t{int(done)}\n")

