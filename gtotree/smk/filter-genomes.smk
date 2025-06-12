import os
import shutil
from gtotree.utils.general import (read_run_data,
                                   write_run_data,
                                   touch)
from gtotree.utils.seqs import filter_seqs_by_genome_ids

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

genome_ids_to_remove = {gd.id for gd in run_data.get_all_input_genomes_due_for_SCG_min_hit_filtering()}

SCG_dict = {SCG.id: SCG for SCG in run_data.get_all_SCG_targets_remaining()}
SCG_ids = list(SCG_dict.keys())

rule all:
    input:
        expand(f"{run_data.found_SCG_seqs_dir}/{{SCG}}.genome-filtered", SCG=SCG_ids)
    run:
        run_data.genomes_filtered_for_min_SCG_hits = True
        for f in input:
            os.remove(f)

        write_run_data(run_data)


rule filter_genomes:
    output:
        f"{run_data.found_SCG_seqs_dir}/{{SCG}}.genome-filtered"
    run:
        inpath = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-gene-filtered{run_data.general_ext}"
        outpath = run_data.found_SCG_seqs_dir + f"/{wildcards.SCG}-genome-filtered{run_data.general_ext}"

        if len(genome_ids_to_remove) == 0:
            shutil.copy(inpath, outpath)
        else:
            filter_seqs_by_genome_ids(inpath, genome_ids_to_remove, outpath)

        touch(output[0])
