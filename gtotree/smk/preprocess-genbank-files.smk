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
            gf = genbank_dict[gb_basename]
            path = gf.full_path
            status_path = f"{run_data.genbank_processing_dir}/{gb_basename}.done"
            with open(status_path, 'r') as f:
                for line in f:
                    gb_file, status, prodigal_used, was_gzipped = line.strip().split('\t')

                    if int(status):
                        if int(was_gzipped):
                            gf.full_path = path[:-3]
                            gf.basename = os.path.basename(gf.full_path)
                            gf.mark_was_gzipped()
                        gf.mark_done()
                    else:
                        gf.mark_removed()

                    if int(prodigal_used):
                        gf.mark_prodigal_used()
                        run_data.tools_used.prodigal_used = True

        write_run_data(run_data)


rule process_genbank_files:
    output:
        f"{run_data.genbank_processing_dir}/{{gb_file}}.done"
    run:
        gf = genbank_dict[wildcards.gb_file]
        path = gf.full_path
        path, was_gzipped = gunzip_if_needed(path)
        if was_gzipped:
            gf.full_path = path
            gf.basename = os.path.basename(path)
        done = extract_filter_and_rename_cds_amino_acids_from_gb(gf.basename, path, run_data)

        if not done:
            extract_fasta_from_gb(gf.basename, path, run_data)
            done = run_prodigal(gf.basename, run_data, path, "genbank")
            prodigal_used = True
            if done:
                filter_and_rename_fasta(gf.basename, run_data, run_data.genbank_processing_dir)
            else:
                prodigal_used = False
        else:
            prodigal_used = False

        with open(output[0], 'w') as f:
            f.write(f'{wildcards.gb_file}\t{int(done)}\t{int(prodigal_used)}\t{int(was_gzipped)}\n')
