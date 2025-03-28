import os
from gtotree.utils.general import (read_run_data,
                                   write_run_data,
                                   touch)
from gtotree.utils.seqs import filter_seqs_by_genome_ids

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

genome_ids_to_remove = {gd.id for gd in run_data.get_all_input_genomes_due_for_SCG_min_hit_filtering()}

rule all:
    input:
        expand(f"{run_data.found_SCG_seqs_dir}/{{SCG}}.genome-filtered", SCG=run_data.remaining_SCG_targets)
    run:
        run_data.genomes_filtered_for_min_SCG_hits = True
        for f in input:
            os.remove(f)

        write_run_data(run_data)


rule filter_genomes:
    output:
        f"{run_data.found_SCG_seqs_dir}/{{SCG}}.genome-filtered"
    run:
        path = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-gene-filtered.fasta"
        filter_seqs_by_genome_ids(path, genome_ids_to_remove)
        touch(output[0])
