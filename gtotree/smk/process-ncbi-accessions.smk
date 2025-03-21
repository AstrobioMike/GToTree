from gtotree.utils.general import (read_run_data,
                                   write_run_data,
                                   read_args,
                                   run_prodigal,
                                   touch)
from gtotree.utils.processing_genomes import prepare_accession

run_data = read_run_data(config['run_data_path'])

accessions = run_data.ncbi_accessions

rule all:
    input:
        expand(f"{run_data.ncbi_downloads_dir}/{{acc}}.done", acc=accessions)
    run:
        for file in input:
            with open(file, 'r') as f:
                for line in f:
                    acc, status = line.strip().split('\t')
                    status = int(status)

                    if status:
                        run_data.add_done_ncbi_accession(acc)
                    else:
                        run_data.remove_ncbi_accession(acc)

        write_run_data(run_data)


rule process_ncbi_accessions:
    output:
        f"{run_data.ncbi_downloads_dir}/{{acc}}.done"
    run:
        done, nt = prepare_accession(wildcards.acc, run_data)
        if done and nt:
            print(f"\n\n DOING {wildcards.acc} \n\n")
            done = run_prodigal(wildcards.acc, run_data, "ncbi")
        # gtt-filter
        # gtt-rename
        # done = filter and rename seqs()
        # more?

        with open(output[0], 'w') as f:
            f.write(f'{wildcards.acc}\t{int(done)}\n')
