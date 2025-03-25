import os
from gtotree.utils.general import (read_run_data,
                                   write_run_data,
                                   read_args,
                                   gunzip_if_needed,
                                   run_prodigal,
                                   touch)
from gtotree.utils.seqs import (filter_and_rename_fasta)

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

fasta_dict = {gf.basename: gf for gf in run_data.fasta_files}
fasta_basenames = list(fasta_dict.keys())

rule all:
    input:
        expand(f"{run_data.fasta_processing_dir}/{{fasta_file}}.done", fasta_file=fasta_basenames)
    run:
        for fasta_basename in fasta_basenames:
            fasta = fasta_dict[fasta_basename]
            path = fasta.full_path
            status_path = f"{run_data.fasta_processing_dir}/{fasta_basename}.done"
            with open(status_path, 'r') as f:
                for line in f:
                    fasta_file, status, was_gzipped = line.strip().split('\t')

                    if int(status):
                        if int(was_gzipped):
                            fasta.full_path = path[:-3]
                            fasta.basename = os.path.basename(fasta.full_path)
                            fasta.mark_was_gzipped()
                        fasta.mark_done()
                    else:
                        fasta.mark_removed()

        write_run_data(run_data)


rule process_fasta_files:
    output:
        f"{run_data.fasta_processing_dir}/{{fasta_file}}.done"
    run:
        fasta = fasta_dict[wildcards.fasta_file]
        path = fasta.full_path
        path, was_gzipped = gunzip_if_needed(path)

        if was_gzipped:
            fasta.full_path = path
            fasta.basename = os.path.basename(path)

        done = run_prodigal(fasta.basename, run_data, fasta.full_path, "fasta")

        if done:
            filter_and_rename_fasta(fasta.basename, run_data, run_data.fasta_processing_dir)

        with open(output[0], 'w') as f:
            f.write(f'{wildcards.fasta_file}\t{int(done)}\t{int(was_gzipped)}\n')
