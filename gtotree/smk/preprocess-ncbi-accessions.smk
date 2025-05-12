from gtotree.utils.general import (read_run_data,
                                   write_run_data)
from gtotree.main_stages.preprocessing_genomes import run_prodigal, prepare_accession
from gtotree.utils.seqs import filter_and_rename_fasta

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

accession_dict = {gd.id: gd for gd in run_data.ncbi_accs if gd.was_found and not gd.preprocessing_done and not gd.removed}
accessions = list(accession_dict.keys())

rule all:
    input:
        expand(f"{run_data.ncbi_downloads_dir}/{{acc}}.done", acc=accessions)
    run:
        for acc in accessions:
            acc_gd = accession_dict[acc]
            status_path = f"{run_data.ncbi_downloads_dir}/{acc}.done"
            with open(status_path, 'r') as f:
                for line in f:
                    acc, status, downloaded, prodigal_used, final_AA_path, num_genes = line.strip().split('\t')

                    acc_gd.num_genes = int(num_genes)

                    if int(status):
                        acc_gd.mark_preprocessing_done()
                        acc_gd.final_AA_path = final_AA_path
                    else:
                        acc_gd.mark_removed()

                    if int(downloaded):
                        acc_gd.was_downloaded = True
                    else:
                        acc_gd.was_downloaded = False
                        acc_gd.reason_removed = "NCBI download failed"

                    acc_gd.prodigal_used = True if int(prodigal_used) else False

        write_run_data(run_data)


rule process_ncbi_accessions:
    output:
        f"{run_data.ncbi_downloads_dir}/{{acc}}.done"
    run:
        acc_gd = accession_dict[wildcards.acc]
        done, nt = prepare_accession(wildcards.acc, run_data)
        downloaded = True if done else False
        if done and nt:
            done = run_prodigal(acc_gd.id, run_data, group = "ncbi")
            prodigal_used = True
        else:
            prodigal_used = False

        if done:
            done, final_AA_path, num_genes = filter_and_rename_fasta(acc_gd.id, run_data, run_data.ncbi_downloads_dir)
        else:
            final_AA_path = None

        with open(output[0], 'w') as f:
            f.write(f'{wildcards.acc}\t{int(done)}\t{int(downloaded)}\t{int(prodigal_used)}\t{final_AA_path}\t{num_genes}\n')
