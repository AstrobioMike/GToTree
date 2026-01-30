import os
import json
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
        expand(f"{run_data.fasta_processing_dir}/{{fasta_file}}.json", fasta_file=fasta_basenames)
    run:
        for fasta_basename in fasta_basenames:
            fasta = fasta_dict[fasta_basename]
            path = fasta.full_path
            status_path = f"{run_data.fasta_processing_dir}/{fasta_basename}.json"

            if not os.path.exists(status_path):
                raise FileNotFoundError(f"Expected status file not found: {status_path}")

            with open(status_path, 'r') as fh:
                data = json.load(fh)

            done = bool(data.get('done', False))
            was_gzipped = bool(data.get('was_gzipped', False))
            prodigal_used = bool(data.get('prodigal_used', False))
            final_AA_path = data.get('final_AA_path')
            final_nt_path = data.get('final_nt_path')
            num_genes = int(data.get('num_genes', 0) or 0)

            fasta.num_genes = num_genes

            if done:
                if was_gzipped:
                    fasta.mark_was_gzipped()
                fasta.mark_preprocessing_done()
                fasta.final_AA_path = final_AA_path
                fasta.final_nt_path = final_nt_path
            else:
                fasta.mark_removed("fasta-file processing failed")
                fasta.preprocessing_failed = True

            if prodigal_used:
                run_data.tools_used.prodigal_used = True

        write_run_data(run_data)


rule process_fasta_files:
    output:
        f"{run_data.fasta_processing_dir}/{{fasta_file}}.json"
    run:
        fasta = fasta_dict[wildcards.fasta_file]
        path, was_gzipped = gunzip_if_needed(fasta.full_path)

        done = run_prodigal(fasta.id, run_data, path, "fasta")

        if was_gzipped:
            os.remove(path)

        if done:
            done, final_AA_path, num_genes, final_nt_path = filter_and_rename_fasta(fasta.id, run_data, run_data.fasta_processing_dir)
        else:
            final_AA_path = None
            final_nt_path = None
            num_genes = 0

        out_obj = {
            "input": wildcards.fasta_file,
            "done": bool(done),
            "was_gzipped": bool(was_gzipped),
            "prodigal_used": bool(done),
            "final_AA_path": final_AA_path,
            "num_genes": int(num_genes or 0),
            "final_nt_path": final_nt_path,
        }

        tmp_path = output[0] + ".tmp"
        with open(tmp_path, 'w') as fh:
            json.dump(out_obj, fh, indent=2, sort_keys=True)
        os.replace(tmp_path, output[0])
