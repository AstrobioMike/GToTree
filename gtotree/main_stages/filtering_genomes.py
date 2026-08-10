import shutil
from gtotree.utils.misc.messaging import (report_message, report_processing_stage,
                                     report_genome_filtering_update)
from gtotree.utils.misc.seqs import check_target_SCGs_have_seqs, filter_seqs_by_genome_ids
from gtotree.utils.misc.general import (write_run_data,
                                   required_count,
                                   run_pooled_stage)

def filter_genomes(args, run_data):

    report_processing_stage("filter-genomes", run_data)
    cutoff = f"{args.genome_hits_cutoff * 100:.0f}"
    if not args.best_hit_mode:
        message = f"Keeping genomes with single hits to at least {cutoff}% of the remaining target-SCGs."
    else:
        message = f"Keeping genomes with hits to at least {cutoff}% of the remaining target-SCGs."
    report_message(message, ii="    ", si="    ", width=80)

    if not run_data.genomes_filtered_for_min_SCG_hits:

        genomes = run_data.get_all_input_genomes_for_filtering()
        num_remaining_SCG_targets = len([SCG.id for SCG in run_data.get_all_SCG_targets_remaining()])
        min_num_SCG_hits = required_count(num_remaining_SCG_targets,
                                          args.genome_hits_cutoff)

        # `or 0` because num_SCG_hits_after_filtering is None until filter_genes has
        # counted; a genome that reached here uncounted has no surviving hits
        genome_ids_to_filter_out = [genome.id for genome in genomes
                                    if (genome.num_SCG_hits_after_filtering or 0) < min_num_SCG_hits]

        for genome in genomes:
            if genome.id in genome_ids_to_filter_out:
                reason = "too few SCG hits" if args.best_hit_mode else "too few unique SCG hits"
                genome.mark_removed(reason)

        genome_ids_to_remove = {gd.id for gd in run_data.get_all_input_genomes_due_for_SCG_min_hit_filtering()}
        scgs = run_data.get_all_SCG_targets_remaining()

        def worker(scg, run_data):
            inpath = run_data.found_SCG_seqs_dir + f"/{scg.id}-gene-filtered{run_data.general_ext}"
            outpath = run_data.found_SCG_seqs_dir + f"/{scg.id}-genome-filtered{run_data.general_ext}"
            # per-SCG file only; no shared state touched here
            if len(genome_ids_to_remove) == 0:
                shutil.copy(inpath, outpath)
            else:
                filter_seqs_by_genome_ids(inpath, genome_ids_to_remove, outpath)
            return

        def apply_result(scg, result, run_data):
            pass

        run_data = run_pooled_stage(scgs, worker, apply_result, args, run_data)

        run_data.genomes_filtered_for_min_SCG_hits = True
        write_run_data(run_data)
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
