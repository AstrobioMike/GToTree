import os
import pandas as pd
from gtotree.utils.general import (read_run_data,
                                   write_run_data)
from gtotree.main_stages.additional_pfam_searching import (run_pfam_search,
                                                           get_pfam_counts,
                                                           write_out_tmp_pfam_hits)

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

genome_dict = {gd.id: gd for gd in run_data.get_all_input_genomes_for_pfam_search()}
genome_ids = list(genome_dict.keys())

rule all:
    input:
        expand(f"{run_data.tmp_pfam_results_dir}/{{ID}}/done", ID=genome_ids)
    run:
        # start_combined_pfam_hit_count_tab(run_data)
        cols = ['assembly_id', 'total_gene_count'] + run_data.found_pfam_targets
        pfam_counts_df = pd.DataFrame(0, index=[], columns=cols)

        for ID in genome_ids:
            genome = genome_dict[ID]
            status_path = f"{run_data.tmp_pfam_results_dir}/{ID}/done"

            with open(status_path, 'r') as f:
                for line in f:
                    (ID, pfam_search_failed) = line.strip().split('\t')

                    if int(pfam_search_failed):
                        genome.mark_pfam_search_failed()
                    else:
                        genome.mark_pfam_search_done()

            # adding info to pfam_counts_df
            pfam_results_txt = f"{run_data.pfam_results_dir}/individual-genome-results/{ID}/pfam-results.txt"
            counts_list = get_pfam_counts(run_data.found_pfam_targets, pfam_results_txt)
            print("Counts list for genome ID:", ID)
            print(counts_list)
            num_genes = genome.num_genes

            pfam_counts_df.loc[len(pfam_counts_df)] = [ID, num_genes] + counts_list

        # writing out the combined pfam counts table
        pfam_counts_df.to_csv(f"{run_data.pfam_results_dir}/pfam-hit-counts.tsv", sep='\t', index=False)

        write_run_data(run_data)


rule search_pfams:
    output:
        f"{run_data.tmp_pfam_results_dir}/{{ID}}/done"
    run:
        genome = genome_dict[wildcards.ID]
        AA_path = genome.final_AA_path

        base_outpath = f"{run_data.pfam_results_dir}/individual-genome-results/{wildcards.ID}"
        pfam_search_failed = run_pfam_search(run_data.all_pfam_targets_hmm_path,
                                         base_outpath,
                                         AA_path,
                                         run_data.num_hmm_cpus)

        if not pfam_search_failed:
            # writing out fastas of hits for each pfam
            results_txt = f"{run_data.pfam_results_dir}/individual-genome-results/{wildcards.ID}/pfam-results.txt"
            out_base = f"{run_data.tmp_pfam_results_dir}/{wildcards.ID}/"
            write_out_tmp_pfam_hits(results_txt, out_base, AA_path)

        with open(f"{output}", 'w') as done_file:
            done_file.write(f"{wildcards.ID}\t{int(pfam_search_failed)}\n")
