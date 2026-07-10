from gtotree.utils.messaging import (report_message, report_processing_stage,
                                     report_genome_filtering_update)
from gtotree.utils.seqs import check_target_SCGs_have_seqs
from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   run_snakemake)

def filter_genomes(args, run_data):

    report_processing_stage("filter-genomes", run_data)
    cutoff = "{:.0f}".format(args.genome_hits_cutoff * 100)
    if not args.best_hit_mode:
        message = f"Keeping those with single hits to at least {cutoff}% of the remaining target-SCGs."
    else:
        message = f"Keeping those with hits to at least {cutoff}% of the remaining target-SCGs."
    report_message(message, ii="    ", si="    ", width=80)

    if not run_data.genomes_filtered_for_min_SCG_hits:

        genomes = run_data.get_all_input_genomes_for_filtering()
        num_remaining_SCG_targets = len([SCG.id for SCG in run_data.get_all_SCG_targets_remaining()])
        # not using round() to avoid banker's rounding, and already checked up front these will always be positive
        min_num_SCG_hits = int(num_remaining_SCG_targets * args.genome_hits_cutoff + 0.5)

        genome_ids_to_filter_out = [genome.id for genome in genomes if genome.num_SCG_hits_after_filtering < min_num_SCG_hits]

        for genome in genomes:
            if genome.id in genome_ids_to_filter_out:
                reason = "too few SCG hits" if args.best_hit_mode else "too few unique SCG hits"
                genome.mark_removed(reason)

        write_run_data(run_data)
        snakefile = get_snakefile_path("filter-genomes.smk")
        description = "Filtering genomes"

        run_snakemake(snakefile, num_remaining_SCG_targets, args, run_data, description)

        run_data = read_run_data(run_data.run_data_path)
        capture_removed_genomes(run_data)

        run_data = check_target_SCGs_have_seqs(run_data, f"-genome-filtered{run_data.general_ext}")

    report_genome_filtering_update(run_data)

    return run_data


def capture_removed_genomes(run_data):

    removed = run_data.get_all_input_genomes_due_for_SCG_min_hit_filtering()

    if len(removed) > 0:

        out_path = run_data.run_files_dir + "/genomes-removed-for-too-few-SCG-hits.tsv"

        with open(out_path, "w") as fail_file:

            fail_file.write("assembly_id\ttotal_SCG_hits\tunique_SCG_hits\tnum_SCG_hits_after_filtering\n")
            for genome in removed:
                total_hits = getattr(genome, 'num_SCG_hits', 0) or 0
                unique_hits = getattr(genome, 'num_unique_SCG_hits', 0) or 0
                hits_after = getattr(genome, 'num_SCG_hits_after_filtering', 0) or 0
                fail_file.write(f"{genome.id}\t{int(total_hits)}\t{int(unique_hits)}\t{int(hits_after)}\n")
