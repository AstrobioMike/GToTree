import os
from gtotree.utils.general import (read_run_data,
                                   write_run_data)
from gtotree.utils.seqs import filter_seqs_by_length

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

genome_dict = {gd.id: gd for gd in run_data.get_all_input_genomes_for_filtering()}
genome_ids = list(genome_dict.keys())

SCG_dict = {SCG.id: SCG for SCG in run_data.get_all_SCG_targets_remaining()}
SCG_ids = list(SCG_dict.keys())

rule all:
    input:
        expand(f"{run_data.found_SCG_seqs_dir}/{{SCG}}.gene-filtered", SCG=SCG_ids)
    run:
        count_dict = {genome_id: 0 for genome_id in genome_ids}
        for SCG in SCG_ids:
            with open(f"{run_data.found_SCG_seqs_dir}/{SCG}.gene-filtered", 'r') as f:
                for line in f:
                    genome_id = line.strip()
                    count_dict[genome_id] += 1

            SCG_target = SCG_dict[SCG]
            SCG_target.gene_length_filtered = True

        for genome_id, count in count_dict.items():
            genome = genome_dict[genome_id]
            genome.num_SCG_hits_after_filtering = count

        for f in input:
            os.remove(f)

        write_run_data(run_data)


rule filter_genes:
    output:
        f"{run_data.found_SCG_seqs_dir}/{{SCG}}.gene-filtered"
    run:
        path = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}.fasta"
        print(f"\n\n\n{path}\n\n\n")
        genomes_with_hits_after_filtering = filter_seqs_by_length(path, run_data.seq_length_cutoff)
        with open(output[0], 'w') as f:
            for genome in genomes_with_hits_after_filtering:
                f.write(f'{genome}\n')
