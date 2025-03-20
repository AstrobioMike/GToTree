from gtotree.utils.general import (read_run_data,
                                   write_run_data,
                                   read_args,
                                   touch)
from gtotree.utils.processing_genomes import prepare_accession

args = read_args(config['args_path'])
run_data = read_run_data(args)

accessions = run_data.ncbi_accessions

rule all:
    input:
        expand(f"{args.ncbi_downloads_dir}/{{acc}}.done", acc=accessions)
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

        write_run_data(run_data, args)


rule process_ncbi_accessions:
    output:
        f"{args.ncbi_downloads_dir}/{{acc}}.done"
    run:
        done = prepare_accession(wildcards.acc, args, run_data)
        # done = filter and rename seqs
        # more?

        with open(output[0], 'w') as f:
            f.write(f'{wildcards.acc}\t{int(done)}\n')
