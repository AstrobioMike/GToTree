from gtotree.utils.general import (read_run_data,
                                   write_run_data,
                                   read_args,
                                   run_prodigal)
from gtotree.utils.seqs import filter_and_rename_fasta
from gtotree.utils.preprocessing_genomes import prepare_accession

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

accessions = run_data.remaining_ncbi_accessions

rule all:
    input:
        expand(f"{run_data.ncbi_downloads_dir}/{{acc}}.done", acc=accessions)
    run:
        for file in input:
            with open(file, 'r') as f:
                for line in f:
                    acc, status, downloaded, prodigal_used = line.strip().split('\t')

                    if int(status):
                        run_data.add_done_ncbi_accession(acc)
                    else:
                        run_data.remove_ncbi_accession(acc)

                    if not int(downloaded):
                        run_data.add_ncbi_acc_not_downloaded(acc)

                    if int(prodigal_used):
                        run_data.tools_used.prodigal_used = True

        write_run_data(run_data)


rule process_ncbi_accessions:
    output:
        f"{run_data.ncbi_downloads_dir}/{{acc}}.done"
    run:
        done, nt = prepare_accession(wildcards.acc, run_data)
        if done and nt:
            done = run_prodigal(wildcards.acc, run_data, "ncbi")
            prodigal_used = True
        else:
            prodigal_used = False

        if done:
            downloaded = True
            filter_and_rename_fasta(wildcards.acc, run_data, run_data.ncbi_downloads_dir)
        else:
            downloaded = False

        with open(output[0], 'w') as f:
            f.write(f'{wildcards.acc}\t{int(done)}\t{int(downloaded)}\t{int(prodigal_used)}\n')
