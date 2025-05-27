import os
import pandas as pd
from gtotree.utils.general import (read_run_data,
                                   write_run_data)
from gtotree.main_stages.additional_ko_searching import (run_ko_search,
                                                         get_ko_counts,
                                                         write_out_tmp_ko_hits)

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

genome_dict = {gd.id: gd for gd in run_data.get_all_input_genomes_for_ko_search()}
genome_ids = list(genome_dict.keys())

rule all:
    input:
        expand(f"{run_data.tmp_ko_results_dir}/{{ID}}/done", ID=genome_ids)
    run:
        # start_combined_ko_hit_count_tab(run_data)
        cols = ['assembly_id', 'total_gene_count'] + run_data.found_ko_targets
        ko_counts_df = pd.DataFrame(0, index=[], columns=cols)

        for ID in genome_ids:
            genome = genome_dict[ID]
            status_path = f"{run_data.tmp_ko_results_dir}/{ID}/done"

            with open(status_path, 'r') as f:
                for line in f:
                    (ID, ko_search_failed) = line.strip().split('\t')

                    if int(ko_search_failed):
                        genome.mark_ko_search_failed()
                    else:
                        genome.mark_ko_search_done()

            # adding info to ko_counts_df
            ko_results_tsv = f"{run_data.ko_results_dir}/individual-genome-results/{ID}/kofamscan-results.tsv"
            counts_list = get_ko_counts(run_data.found_ko_targets, ko_results_tsv)
            num_genes = genome.num_genes

            ko_counts_df.loc[len(ko_counts_df)] = [ID, num_genes] + counts_list

        # writing out the combined KO counts table
        ko_counts_df.to_csv(f"{run_data.ko_results_dir}/ko-hit-counts.tsv", sep='\t', index=False)

        write_run_data(run_data)


rule search_kos:
    output:
        f"{run_data.tmp_ko_results_dir}/{{ID}}/done"
    run:
        genome = genome_dict[wildcards.ID]
        AA_path = genome.final_AA_path

        base_out_path = f"{run_data.ko_results_dir}/individual-genome-results/{wildcards.ID}/"
        ko_search_failed = run_ko_search(run_data.target_ko_profiles_dir,
                                         run_data.target_kos_tsv,
                                         base_out_path,
                                         AA_path)

        if not ko_search_failed:
            # writing out fastas of hits for each KO
            results_tsv = f"{run_data.ko_results_dir}/individual-genome-results/{wildcards.ID}/kofamscan-results.tsv"
            out_base = f"{run_data.tmp_ko_results_dir}/{wildcards.ID}/"
            write_out_tmp_ko_hits(results_tsv, out_base, AA_path)

        with open(f"{output}", 'w') as done_file:
            done_file.write(f"{wildcards.ID}\t{int(ko_search_failed)}\n")
