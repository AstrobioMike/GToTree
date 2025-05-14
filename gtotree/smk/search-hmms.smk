import os
from gtotree.utils.general import (read_run_data,
                                   write_run_data)
from gtotree.main_stages.hmm_searching import (run_hmm_search,
                                          parse_hmmer_results,
                                          get_seqs,
                                          write_out_SCG_hit_seqs,
                                          start_combined_SCG_hit_count_tab,
                                          add_to_combined_SCG_hit_count_tab)

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

genome_dict = {gd.id: gd for gd in run_data.get_all_input_genomes_for_hmm_search()}
genome_ids = list(genome_dict.keys())

rule all:
    input:
        expand(f"{run_data.hmm_results_dir}/{{ID}}/done", ID=genome_ids)
    run:
        start_combined_SCG_hit_count_tab(run_data)
        for ID in genome_ids:
            genome = genome_dict[ID]
            status_path = f"{run_data.hmm_results_dir}/{ID}/done"

            with open(status_path, 'r') as f:
                for line in f:
                    (ID, hmm_search_failed, extract_seqs_failed,
                    num_SCG_hits, num_unique_SCG_hits) = line.strip().split('\t')

                    if int(hmm_search_failed):
                        genome.mark_hmm_search_failed()
                        genome.mark_removed("HMM search failed")
                        genome.num_SCG_hits = 0
                    else:
                        add_to_combined_SCG_hit_count_tab(genome.id, run_data)
                        if int(extract_seqs_failed):
                            genome.mark_extract_seqs_failed()
                            genome.num_SCG_hits = 0
                            genome.mark_removed("extracting sequences after HMM search failed")
                        else:
                            write_out_SCG_hit_seqs(genome.id, run_data)
                            genome.mark_hmm_search_done()
                            genome.num_SCG_hits = int(num_SCG_hits)
                            genome.num_unique_SCG_hits = int(num_unique_SCG_hits)

        write_run_data(run_data)


rule search_hmms:
    output:
        f"{run_data.hmm_results_dir}/{{ID}}/done"
    run:
        genome = genome_dict[wildcards.ID]
        AA_path = genome.final_AA_path

        out_dir = f"{run_data.hmm_results_dir}/{wildcards.ID}"
        hmm_out_path = f"{out_dir}/SCG-hits-hmm.txt"
        hmm_search_failed = run_hmm_search(wildcards.ID, run_data, AA_path, hmm_out_path)

        if not hmm_search_failed:

            (dict_of_hit_counts, dict_of_hit_gene_ids,
            num_SCG_hits, num_unique_SCG_hits) = parse_hmmer_results(f"{out_dir}/SCG-hits-hmm.txt",
                                                                                         run_data)

            with open(f"{out_dir}/SCG-hit-counts.txt", 'w') as f:
                f.write(f"{genome.id}")
                for count in dict_of_hit_counts.values():
                    f.write(f"\t{count}")
                f.write("\n")

            hit_seqs_dict, extract_seqs_failed = get_seqs(dict_of_hit_gene_ids, AA_path)

            if not extract_seqs_failed:
                with open(f"{out_dir}/SCG-hits.fasta", 'w') as f:
                    for gene_id, seq in hit_seqs_dict.items():
                        if seq is not None:
                            f.write(f">{gene_id}\n{seq}\n")

        with open(output[0], 'w') as f:
            f.write(f"{wildcards.ID}\t{int(hmm_search_failed)}\t{int(extract_seqs_failed)}\t{int(num_SCG_hits)}\t{int(num_unique_SCG_hits)}\n")
