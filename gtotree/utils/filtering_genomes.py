from gtotree.utils.messaging import (report_processing_stage,
                                     report_genome_filtering_update)
from gtotree.utils.seqs import check_target_SCGs_have_seqs
from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   run_snakemake)

def filter_genomes(args, run_data):

    report_processing_stage("filter-genomes")
    cutoff = "{:.0f}".format(args.genome_hits_cutoff * 100)
    if not args.best_hit_mode:
        print(f"\n    Keeping those with single hits to at least {cutoff}% of the total targeted SCGs.")
    else:
        print(f"\n       Keeping those with hits to at least {cutoff}% of the total targeted SCGs.")

    genomes = run_data.get_all_input_genomes_for_filtering()
    num_remaining_SCG_targets = len([SCG.id for SCG in run_data.get_all_SCG_targets_remaining()])
    min_num_SCG_hits = round(num_remaining_SCG_targets * args.genome_hits_cutoff)
    genome_ids_to_filter_out = [genome.id for genome in genomes if genome.num_SCG_hits < min_num_SCG_hits]

    for genome in genomes:
        if genome.id in genome_ids_to_filter_out:
            genome.removed = True
            genome.reason_removed = "too few SCG hits"

    write_run_data(run_data)
    snakefile = get_snakefile_path("filter-genomes.smk")
    description = "Filtering genomes"

    run_snakemake(snakefile, num_remaining_SCG_targets, args, run_data, description)

    run_data = read_run_data(run_data.run_data_path)
    capture_removed_genomes(run_data)

    run_data = check_target_SCGs_have_seqs(run_data, "-genome-filtered.fasta")

    report_genome_filtering_update(run_data)

    return run_data


def capture_removed_genomes(run_data):
    if len(run_data.get_all_input_genomes_due_for_SCG_min_hit_filtering()) > 0:
        with open(run_data.run_files_dir + "/genomes-removed-for-too-few-SCG-hits.txt", "w") as fail_file:
            for genome in run_data.get_all_input_genomes_due_for_SCG_min_hit_filtering():
                fail_file.write(genome.id + "\n")
