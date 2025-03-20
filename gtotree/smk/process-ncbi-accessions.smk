from gtotree.utils.general import (read_genome_data,
                                   write_genome_data,
                                   read_args,
                                   touch)
from gtotree.utils.processing_genomes import prepare_accession

genome_data = read_genome_data(config['genome_data_path'])
args = read_args(config['args_path'])

accessions = genome_data.ncbi_accessions

rule all:
    input:
        expand(f"{args.ncbi_downloads_dir}/{{acc}}.done", acc=accessions)
    run:
        write_genome_data(genome_data, args)


rule process_ncbi_accessions:
    output:
        f"{args.ncbi_downloads_dir}/{{acc}}.done"
    run:
        downloaded = prepare_accession(wildcards.acc, args, genome_data)
        if not downloaded:
            genome_data.remove_ncbi_accession(wildcards.acc)
        touch(output[0])
