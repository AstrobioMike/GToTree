import os
from gtotree.utils.general import (read_run_data,
                                   write_run_data,
                                   gunzip_if_needed)
from gtotree.main_stages.preprocessing_genomes import run_prodigal
from gtotree.utils.seqs import filter_and_rename_fasta

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

fasta_dict = {gd.basename: gd for gd in run_data.fasta_files if not gd.preprocessing_done and not gd.removed}
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
                    fasta_file, status, was_gzipped, final_AA_path, num_genes = line.strip().split('\t')

                    fasta.num_genes = int(num_genes)

                    if int(status):
                        if int(was_gzipped):
                            fasta.mark_was_gzipped()
                        fasta.mark_preprocessing_done()
                        fasta.final_AA_path = final_AA_path
                    else:
                        fasta.mark_removed("fasta-file processing failed")
                        fasta.preprocessing_failed = True


        run_data.tools_used.prodigal_used = True
        write_run_data(run_data)


rule process_fasta_files:
    output:
        f"{run_data.fasta_processing_dir}/{{fasta_file}}.done"
    run:
        fasta = fasta_dict[wildcards.fasta_file]
        path, was_gzipped = gunzip_if_needed(fasta.full_path)

        done = run_prodigal(fasta.id, run_data, path, "fasta")

        if was_gzipped:
            os.remove(path)

        if done:
            done, final_AA_path, num_genes = filter_and_rename_fasta(fasta.id, run_data, run_data.fasta_processing_dir)
        else:
            final_AA_path = None

        with open(output[0], 'w') as f:
            f.write(f'{wildcards.fasta_file}\t{int(done)}\t{int(was_gzipped)}\t{final_AA_path}\t{num_genes}\n')
