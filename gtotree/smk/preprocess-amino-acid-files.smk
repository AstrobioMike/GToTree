import os
import json
from gtotree.utils.general import (read_run_data,
                                   write_run_data,
                                   gunzip_if_needed)
from gtotree.utils.seqs import filter_and_rename_fasta

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

AA_dict = {gd.basename: gd for gd in run_data.amino_acid_files if not gd.preprocessing_done and not gd.removed}
AA_basenames = list(AA_dict.keys())

rule all:
    input:
        expand(f"{run_data.AA_processing_dir}/{{AA_file}}.json", AA_file=AA_basenames)
    run:
        for AA_basename in AA_basenames:
            AA = AA_dict[AA_basename]
            path = AA.full_path
            status_path = f"{run_data.AA_processing_dir}/{AA_basename}.json"

            if not os.path.exists(status_path):
                raise FileNotFoundError(f"Expected status file not found: {status_path}")

            with open(status_path, 'r') as fh:
                data = json.load(fh)

            done = bool(data.get('done', False))
            was_gzipped = bool(data.get('was_gzipped', False))
            final_AA_path = data.get('final_AA_path')
            final_nt_path = data.get('final_nt_path')
            num_genes = int(data.get('num_genes', 0) or 0)

            AA.num_genes = num_genes

            if done:
                if was_gzipped:
                    AA.mark_was_gzipped()
                AA.mark_preprocessing_done()
                AA.final_AA_path = final_AA_path
                AA.final_nt_path = final_nt_path
            else:
                AA.mark_removed("amino-acid-file processing failed")
                AA.preprocessing_failed = True

        write_run_data(run_data)


rule process_AA_files:
    output:
        f"{run_data.AA_processing_dir}/{{AA_file}}.json"
    run:
        AA = AA_dict[wildcards.AA_file]
        path, was_gzipped = gunzip_if_needed(AA.full_path)

        done, final_AA_path, num_genes, final_nt_path = filter_and_rename_fasta(AA.id, run_data, path, full_path = True)
        if was_gzipped:
            os.remove(path)

        out_obj = {
            "input": wildcards.AA_file,
            "done": bool(done),
            "was_gzipped": bool(was_gzipped),
            "prodigal_used": False,
            "final_AA_path": final_AA_path,
            "num_genes": int(num_genes or 0),
            "final_nt_path": final_nt_path,
        }

        tmp_path = output[0] + ".tmp"
        with open(tmp_path, 'w') as fh:
            json.dump(out_obj, fh, indent=2, sort_keys=True)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp_path, output[0])
