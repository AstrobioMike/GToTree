import json
import os
from gtotree.utils.general import (read_run_data,
                                   write_run_data)
from gtotree.main_stages.preprocessing_genomes import run_prodigal, prepare_accession
from gtotree.utils.seqs import filter_and_rename_fasta


run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

accession_dict = {gd.id: gd for gd in run_data.ncbi_accs if gd.acc_was_found and not gd.preprocessing_done and not gd.removed}
accessions = list(accession_dict.keys())

rule all:
    input:
        expand(f"{run_data.ncbi_processing_dir}/{{acc}}.json", acc=accessions)
    run:
        for acc in accessions:
            acc_gd = accession_dict[acc]
            status_path = f"{run_data.ncbi_processing_dir}/{acc}.json"

            if not os.path.exists(status_path):
                raise FileNotFoundError(f"Expected status file not found: {status_path}")

            with open(status_path, 'r') as fh:
                data = json.load(fh)

            # normalize/ensure fields
            done = bool(data.get('done') or False)
            downloaded = bool(data.get('downloaded', False))
            prodigal_used = bool(data.get('prodigal_used', False))
            final_AA_path = data.get('final_AA_path')
            final_nt_path = data.get('final_nt_path')
            num_genes = int(data.get('num_genes', 0) or 0)

            acc_gd.num_genes = num_genes

            if done:
                acc_gd.mark_preprocessing_done()
                acc_gd.final_AA_path = final_AA_path
                acc_gd.final_nt_path = final_nt_path
                acc_gd.acc_was_downloaded = downloaded
            else:
                if downloaded:
                    acc_gd.acc_was_downloaded = True
                    acc_gd.mark_removed("acc processing failed after download")
                else:
                    acc_gd.acc_was_downloaded = False
                    acc_gd.mark_removed("acc download failed")

            acc_gd.prodigal_used = prodigal_used

        write_run_data(run_data)


rule process_ncbi_accessions:
    output:
        f"{run_data.ncbi_processing_dir}/{{acc}}.json"
    run:
        acc_gd = accession_dict[wildcards.acc]
        done, nt = prepare_accession(wildcards.acc, run_data)
        downloaded = bool(done)
        if done and nt:
            done = run_prodigal(acc_gd.id, run_data, group = "ncbi")
            prodigal_used = True
        else:
            prodigal_used = False

        if done:
            done, final_AA_path, num_genes, final_nt_path = filter_and_rename_fasta(acc_gd.id, run_data, run_data.ncbi_processing_dir)
        else:
            final_AA_path = None
            final_nt_path = None
            num_genes = 0

        out_obj = {
            "input": wildcards.acc,
            "done": bool(done),
            "downloaded": downloaded,
            "prodigal_used": prodigal_used,
            "final_AA_path": final_AA_path,
            "num_genes": int(num_genes or 0),
            "final_nt_path": final_nt_path,
        }

        tmp_path = output[0] + ".json.tmp"
        with open(tmp_path, 'w') as fh:
            json.dump(out_obj, fh, indent=2, sort_keys=True)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp_path, output[0])
