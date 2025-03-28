import os
from gtotree.utils.general import (read_run_data,
                                   write_run_data)
from gtotree.utils.seqs import filter_seqs_by_length

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

genome_dict = {gd.id: gd for gd in run_data.get_all_input_genomes_for_filtering()}
genome_ids = list(genome_dict.keys())

rule all:
    input:
        expand(f"{run_data.found_SCG_seqs_dir}/{{SCG}}.gene-filtered", SCG=run_data.remaining_SCG_targets)
    run:
        run_data.SCG_hits_filtered = True

        count_dict = {genome_id: 0 for genome_id in genome_ids}
        for SCG in run_data.remaining_SCG_targets:
            with open(f"{run_data.found_SCG_seqs_dir}/{SCG}.gene-filtered", 'r') as f:
                for line in f:
                    genome_id = line.strip()
                    count_dict[genome_id] += 1

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
        genomes_with_hits_after_filtering = filter_seqs_by_length(path, run_data.seq_length_cutoff)
        with open(output[0], 'w') as f:
            for genome in genomes_with_hits_after_filtering:
                f.write(f'{genome}\n')
