import os
from gtotree.utils.general import (read_run_data,
                                   write_run_data,
                                   read_args,
                                   gunzip_if_needed,
                                   run_prodigal,
                                   touch)
from gtotree.utils.seqs import (extract_filter_and_rename_cds_amino_acids_from_gb,
                                extract_fasta_from_gb,
                                filter_and_rename_fasta)

run_data = read_run_data(config['run_data_path'])
if run_data is None:
    raise ValueError("Run data not found")

genbank_dict = {gf.basename: gf for gf in run_data.genbank_files}
genbank_basenames = list(genbank_dict.keys())

rule all:
    input:
        expand(f"{run_data.genbank_processing_dir}/{{gb_file}}.done", gb_file=genbank_basenames)
    run:
        for gb_basename in genbank_basenames:
            gb = genbank_dict[gb_basename]
            path = gb.full_path
            status_path = f"{run_data.genbank_processing_dir}/{gb_basename}.done"
            with open(status_path, 'r') as f:
                for line in f:
                    gb_file, status, prodigal_used, was_gzipped, final_AA_path = line.strip().split('\t')

                    if int(status):
                        if int(was_gzipped):
                            gb.mark_was_gzipped()
                        gb.mark_done()
                        gb.final_AA_path = final_AA_path
                    else:
                        gb.mark_removed()

                    if int(prodigal_used):
                        gb.mark_prodigal_used()
                        run_data.tools_used.prodigal_used = True

        write_run_data(run_data)


rule process_genbank_files:
    output:
        f"{run_data.genbank_processing_dir}/{{gb_file}}.done"
    run:
        gb = genbank_dict[wildcards.gb_file]
        path, was_gzipped = gunzip_if_needed(gb.full_path)

        done, final_AA_path = extract_filter_and_rename_cds_amino_acids_from_gb(gb.id, path, run_data)

        if not done:
            extract_fasta_from_gb(gb.id, path, run_data)
            done = run_prodigal(gb.id, run_data, path, "genbank")
            prodigal_used = True
            if done:
                done, final_AA_path = filter_and_rename_fasta(gb.id, run_data, run_data.genbank_processing_dir)
            else:
                prodigal_used = False
        else:
            prodigal_used = False

        if was_gzipped:
            os.remove(path)

        with open(output[0], 'w') as f:
            f.write(f'{wildcards.gb_file}\t{int(done)}\t{int(prodigal_used)}\t{int(was_gzipped)}\t{final_AA_path}\n')
