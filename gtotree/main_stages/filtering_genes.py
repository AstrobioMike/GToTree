from gtotree.utils.general import (write_run_data,
                                   run_pooled_stage)
from gtotree.utils.messaging import (report_message,
                                     report_processing_stage,
                                     report_SCG_set_filtering_update)
from gtotree.utils.seqs import check_target_SCGs_have_seqs, filter_seqs_by_length

def filter_genes(args, run_data):

    report_processing_stage("filter-genes", run_data)
    cutoff = "{:.0f}".format(run_data.seq_length_cutoff * 100)
    in_genomes_cutoff_for_report = "{:.0f}".format(args.gene_representation_cutoff * 100)

    message = (f"Keeping genes with lengths within {cutoff}% of the median for each gene set, "
               f"and keeping gene sets with hits in at least {in_genomes_cutoff_for_report}% of "
               f"the currently retained genomes.")
    report_message(message, ii="    ", si="    ", width=80)

    scgs_to_filter = run_data.get_all_SCG_targets_remaining_but_not_filtered()

    if len(scgs_to_filter) > 0:
        genome_dict = {gd.id: gd for gd in run_data.get_all_input_genomes_for_filtering()}

        # accumulate per-genome surviving-hit counts across SCGs on the main thread
        count_dict = {genome_id: 0 for genome_id in genome_dict}

        def worker(scg, run_data):
            path = run_data.found_SCG_seqs_dir + f"/{scg.id}{run_data.general_ext}"
            # returns the ids of genomes whose hit survived length-filtering; the
            # length-filtered FASTA is written alongside as a side effect (per-SCG
            # file only -- no shared state touched here)
            return filter_seqs_by_length(path, run_data.seq_length_cutoff)

        def apply_result(scg, genomes_with_hits, run_data):
            scg.gene_length_filtered = True
            scg.num_genomes_with_hits_after_len_filtering = len(genomes_with_hits)
            for genome_id in genomes_with_hits:
                if genome_id in count_dict:
                    count_dict[genome_id] += 1

        run_data = run_pooled_stage(scgs_to_filter, worker, apply_result, args, run_data)

        for genome_id, count in count_dict.items():
            genome_dict[genome_id].num_SCG_hits_after_filtering = count

        write_run_data(run_data)

    run_data = check_target_SCGs_have_seqs(run_data, f"-gene-filtered{run_data.general_ext}")

    total_genomes_remaining = len(run_data.get_all_input_genomes_for_filtering())
    # not using round() to avoid banker's rounding, and already checked up front these will always be positive
    min_genomes_required = int(total_genomes_remaining * args.gene_representation_cutoff + 0.5)

    removed_any = False
    for scg in run_data.get_all_SCG_targets_remaining():
        count = getattr(scg, 'num_genomes_with_hits_after_len_filtering', 0)
        if count < min_genomes_required:
            scg.mark_removed(f"too few genomes with hits")
            removed_any = True

    if removed_any:
        write_run_data(run_data)

    write_out_removed_SCG_targets(run_data)

    report_SCG_set_filtering_update(run_data)

    return run_data


def write_out_removed_SCG_targets(run_data):

    removed_scg_objs = [scg for scg in run_data.SCG_targets if getattr(scg, 'removed', False)]
    if len(removed_scg_objs) > 0:
        out_path = run_data.run_files_dir + "/target-SCGs-dropped-from-analysis.tsv"
        with open(out_path, "w") as fail_file:
            # header
            fail_file.write("target_SCG\treason_removed\n")
            for scg in removed_scg_objs:
                reason = scg.reason_removed or ""
                fail_file.write(f"{scg.id}\t{reason}\n")
