import os
from gtotree.utils.general import (read_run_data,
                                   write_run_data,
                                   read_args,
                                   gunzip_if_needed,
                                   run_prodigal,
                                   touch)
from gtotree.utils.seqs import filter_and_rename_fasta

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

genbank_dict = {os.path.basename(f): f for f in run_data.genbank_files_full_paths}
genbank_basenames = list(genbank_dict.keys())

rule all:
    input:
        expand(f"{run_data.genbank_processing_dir}/{{gb_file}}.done", gb_file=genbank_basenames)
    run:
        for gb_basename in genbank_basenames:
            path = genbank_dict[gb_basename]
            status_path = f"{run_data.genbank_processing_dir}/{gb_basename}.done"
            with open(status_path, 'r') as f:
                for line in f:
                    gb_basename, status, prodigal_used = line.strip().split('\t')

                    if int(status):
                        if gb_basename.endswith('.gz'):
                            run_data.replace_in_list('genbank_files_full_paths', path, path[:-3])
                            path = path[:-3]
                        run_data.add_done_genbank_file(path)
                    else:
                        run_data.remove_genbank_file(path)

                    if int(prodigal_used):
                        run_data.tools_used.prodigal_used = True

        write_run_data(run_data)


rule process_genbank_files:
    output:
        f"{run_data.genbank_processing_dir}/{{gb_file}}.done"
    run:
        path = genbank_dict[wildcards.gb_file]
        path = gunzip_if_needed(path)
        done = True
        prodigal_used = False
        # try to get amino acids
        # if not, get fasta and run prodigal

        #     done = run_prodigal(wildcards.gb_file, run_data, "genbank")
        #     prodigal_used = True
        # else:
        #     prodigal_used = False

        # if done:
        #     downloaded = True
        #     filter_and_rename_fasta(wildcards.acc, run_data)
        # else:
        #     downloaded = False

        with open(output[0], 'w') as f:
            f.write(f'{wildcards.gb_file}\t{int(done)}\t{int(prodigal_used)}\n')
