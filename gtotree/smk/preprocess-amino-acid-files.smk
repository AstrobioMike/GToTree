import os
from gtotree.utils.general import (read_run_data,
                                   write_run_data,
                                   gunzip_if_needed)
from gtotree.utils.seqs import filter_and_rename_fasta

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

AA_dict = {gd.basename: gd for gd in run_data.amino_acid_files}
AA_basenames = list(AA_dict.keys())

rule all:
    input:
        expand(f"{run_data.AA_processing_dir}/{{AA_file}}.done", AA_file=AA_basenames)
    run:
        for AA_basename in AA_basenames:
            AA = AA_dict[AA_basename]
            path = AA.full_path
            status_path = f"{run_data.AA_processing_dir}/{AA_basename}.done"
            with open(status_path, 'r') as f:
                for line in f:
                    AA_file, status, was_gzipped, final_AA_path = line.strip().split('\t')

                    if int(status):
                        if int(was_gzipped):
                            AA.mark_was_gzipped()
                        AA.mark_done()
                        AA.final_AA_path = final_AA_path
                    else:
                        AA.mark_removed()

        write_run_data(run_data)


rule process_AA_files:
    output:
        f"{run_data.AA_processing_dir}/{{AA_file}}.done"
    run:
        AA = AA_dict[wildcards.AA_file]
        path, was_gzipped = gunzip_if_needed(AA.full_path)

        done, final_AA_path = filter_and_rename_fasta(AA.id, run_data, path, full_path = True)
        if was_gzipped:
            os.remove(path)

        with open(output[0], 'w') as f:
            f.write(f'{wildcards.AA_file}\t{int(done)}\t{int(was_gzipped)}\t{final_AA_path}\n')
