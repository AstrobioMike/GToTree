import os
import json
from gtotree.utils.general import (read_run_data,
                                   write_run_data,
                                   gunzip_if_needed)
from gtotree.utils.seqs import (extract_filter_and_rename_cds_amino_acids_from_gb,
                                extract_fasta_from_gb,
                                filter_and_rename_fasta)
from gtotree.main_stages.preprocessing_genomes import run_prodigal

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

genbank_dict = {gd.basename: gd for gd in run_data.genbank_files if not gd.preprocessing_done and not gd.removed}
genbank_basenames = list(genbank_dict.keys())

rule all:
    input:
        expand(f"{run_data.genbank_processing_dir}/{{gb_file}}.json", gb_file=genbank_basenames)
    run:
        for gb_basename in genbank_basenames:
            gb = genbank_dict[gb_basename]
            path = gb.full_path
            status_path = f"{run_data.genbank_processing_dir}/{gb_basename}.json"

            if not os.path.exists(status_path):
                raise FileNotFoundError(f"Expected status file not found: {status_path}")

            with open(status_path, 'r') as fh:
                data = json.load(fh)

            done = bool(data.get('done', False))
            prodigal_used = bool(data.get('prodigal_used', False))
            was_gzipped = bool(data.get('was_gzipped', False))
            final_AA_path = data.get('final_AA_path')
            final_nt_path = data.get('final_nt_path')
            num_genes = int(data.get('num_genes', 0) or 0)

            gb.num_genes = num_genes

            if done:
                if was_gzipped:
                    gb.mark_was_gzipped()
                gb.mark_preprocessing_done()
                gb.final_AA_path = final_AA_path
                gb.final_nt_path = final_nt_path
            else:
                gb.mark_removed("genbank-file processing failed")
                gb.preprocessing_failed = True

            if prodigal_used:
                gb.mark_prodigal_used()
                run_data.tools_used.prodigal_used = True

        write_run_data(run_data)


rule process_genbank_files:
    output:
        f"{run_data.genbank_processing_dir}/{{gb_file}}.json"
    run:
        gb = genbank_dict[wildcards.gb_file]
        path, was_gzipped = gunzip_if_needed(gb.full_path)

        prodigal_used = False
        final_nt_path = None
        final_AA_path = None
        num_genes = 0

        done, final_AA_path, num_genes = extract_filter_and_rename_cds_amino_acids_from_gb(gb.id, path, run_data)

        if not done:
            extract_fasta_from_gb(gb.id, path, run_data)
            done = run_prodigal(gb.id, run_data, path, "genbank")
            prodigal_used = True
            if done:
                done, final_AA_path, num_genes, final_nt_path = filter_and_rename_fasta(gb.id, run_data, run_data.genbank_processing_dir)
            else:
                prodigal_used = False
        else:
            prodigal_used = False

        if was_gzipped:
            os.remove(path)

        out_obj = {
            "input": wildcards.gb_file,
            "done": bool(done),
            "prodigal_used": prodigal_used,
            "was_gzipped": bool(was_gzipped),
            "final_AA_path": final_AA_path,
            "num_genes": int(num_genes or 0),
            "final_nt_path": final_nt_path,
        }

        tmp_path = output[0] + ".tmp"
        with open(tmp_path, 'w') as fh:
            json.dump(out_obj, fh, indent=2, sort_keys=True)
        os.replace(tmp_path, output[0])
