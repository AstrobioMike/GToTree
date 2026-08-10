import tempfile

from gtotree.utils.misc.messaging import (report_processing_stage,
                                     report_ncbi_update,
                                     report_genbank_update,
                                     report_fasta_update,
                                     report_AA_update,
                                     report_genome_processing_update,
                                     report_pfam_searching_update,
                                     report_ko_searching_update)
from gtotree.utils.misc import phase_stats
from gtotree.utils.misc.general import write_run_data, run_pooled_stage, remove_file_if_exists
from gtotree.utils.misc.seqs import check_target_SCGs_have_seqs
from gtotree.utils.misc.stages import PipelineStage, SCGRemovalStage
from gtotree.utils.pfam.pfam_handling import get_additional_pfam_targets
from gtotree.utils.ko.ko_handling import parse_kofamscan_targets
from gtotree.utils.hmms.hmm_searching_engine import press_profiles

from gtotree.utils.misc.processing_genomes import (
    build_base_link_map,
    capture_ncbi_failed_downloads,
    capture_failed_genbank_files,
    capture_failed_fasta_files,
    capture_failed_amino_acid_files,
    _process_one_ncbi_accession,
    _apply_ncbi_accession_status,
    _process_one_genbank_file,
    _apply_genbank_status,
    _process_one_fasta_file,
    _apply_fasta_status,
    _process_one_amino_acid_file,
    _apply_amino_acid_status,
)
from gtotree.utils.hmms.hmm_searching import (
    _hmm_search_worker,
    _apply_hmm_search_result,
    capture_hmm_search_failures,
    rebuild_combined_SCG_outputs,
)
from gtotree.utils.pfam.additional_pfam_searching import (
    _pfam_search_worker,
    combine_all_pfam_hits,
    write_pfam_counts_table,
)
from gtotree.utils.ko.additional_ko_searching import (
    _ko_search_worker,
    combine_all_ko_hits,
    write_ko_counts_table,
    write_out_failed_ko_targets,
)


class SearchPlan:
    """
    Which per-genome searches this run performs, resolved once up front

    The target sets (which Pfams, which KOs) are expensive to work out and identical
    for every genome, so they're resolved here rather than inside the worker.

    `do_scg` is True for the main GToTree run, which always searches the SCG set. The
    `gtt search-pfams` / `gtt search-kos` subcommands reuse this same fused
    preprocess-then-search machinery but have no SCG set and no tree, so they set it
    False
    """

    __slots__ = ("do_pfam", "do_ko", "do_scg", "keep_genome_files",
                 "pressed_scg_base", "pressed_pfam_base")

    def __init__(self, do_pfam, do_ko, keep_genome_files, do_scg=True):
        self.do_pfam = do_pfam
        self.do_ko = do_ko
        self.do_scg = do_scg
        self.keep_genome_files = keep_genome_files
        self.pressed_scg_base = None
        self.pressed_pfam_base = None


def build_search_plan(args, run_data):
    """
    Resolve the target sets and return the plan. Must run before any worker starts
    """
    do_pfam = False
    do_ko = False

    if run_data.target_pfams_file and not run_data.additional_pfam_searching_done:
        run_data = get_additional_pfam_targets(run_data)
        do_pfam = len(run_data.found_pfam_targets) > 0
    elif run_data.target_pfams_file:
        do_pfam = len(run_data.found_pfam_targets) > 0

    if run_data.target_kos_file and not run_data.additional_ko_searching_done:
        run_data = parse_kofamscan_targets(run_data)
        write_out_failed_ko_targets(run_data)
        do_ko = len(run_data.found_ko_targets) > 0
    elif run_data.target_kos_file:
        do_ko = len(run_data.found_ko_targets) > 0

    plan = SearchPlan(do_pfam=do_pfam, do_ko=do_ko, keep_genome_files=bool(args.debug))
    return run_data, plan


def genome_is_fully_processed(gd, plan):
    """
    True when nothing more needs doing for this genome
    """
    if not gd.processing_done:
        return False
    if plan.do_scg and not gd.hmm_search_done:
        return False
    if plan.do_pfam and not gd.pfam_search_done:
        return False
    if plan.do_ko and not gd.ko_search_done:
        return False
    return True


def genomes_needing_processing(source_list, plan):
    return [gd for gd in source_list
            if not gd.removed and not genome_is_fully_processed(gd, plan)]


def _drop_genome_files(status, run_data):
    """
    Delete this genome's FASTAs now that every search that reads them has finished
    """
    for key in ("final_AA_path", "final_nt_path"):
        path = status.get(key)
        if path:
            remove_file_if_exists(path)


def _run_searches(genome, run_data, plan, status):
    """
    Run every configured search for one genome
    """
    aa_path = status.get("final_AA_path")
    nt_path = status.get("final_nt_path")

    results = {}

    if plan.do_pfam:
        results["pfam"] = _pfam_search_worker(genome, run_data, aa_path=aa_path,
                                              pressed_base=plan.pressed_pfam_base)

    if plan.do_ko:
        results["ko"] = _ko_search_worker(genome, run_data, aa_path=aa_path)

    if plan.do_scg:
        results["hmm"] = _hmm_search_worker(genome, run_data,
                                            aa_path=aa_path, nt_path=nt_path,
                                            pressed_base=plan.pressed_scg_base)
    return results


def _apply_searches(genome, searches, run_data, plan):
    """
    Fold the search results into the GenomeData
    """
    if not searches:
        return

    if plan.do_pfam and "pfam" in searches:
        if searches["pfam"].get("pfam_search_failed"):
            genome.mark_pfam_search_failed()
        else:
            genome.mark_pfam_search_done()

    if plan.do_ko and "ko" in searches:
        if searches["ko"].get("ko_search_failed"):
            genome.mark_ko_search_failed()
        else:
            genome.mark_ko_search_done()

    if plan.do_scg and "hmm" in searches:
        _apply_hmm_search_result(genome, searches["hmm"], run_data)


def _fused(preprocess_worker, preprocess_apply, plan):
    """
    Wrap a source-specific processing worker/apply pair into a fused pair that also
    searches the genome and then drops its FASTAs
    """
    def worker(gd, run_data):
        try:
            status = preprocess_worker(gd, run_data)
            if status.get("done"):
                status["searches"] = _run_searches(gd, run_data, plan, status)
                if not plan.keep_genome_files:
                    _drop_genome_files(status, run_data)
            return status
        except BaseException as e:
            # belt and braces: the processing workers already swallow everything,
            # but a raise here would abort the whole stage with some genomes applied
            # and some not
            return {"done": False, "num_genes": 0, "prodigal_used": False,
                    "was_gzipped": False, "final_AA_path": None, "final_nt_path": None,
                    "error": f"{type(e).__name__}: {e}"}

    def apply_result(gd, status, run_data):
        preprocess_apply(gd, status, run_data)
        if status.get("done"):
            _apply_searches(gd, status.get("searches"), run_data, plan)

    return worker, apply_result


def process_genomes(args, run_data):
    """
    Download/preprocess and search every input genome, one genome per worker call
    """
    run_data, plan = build_search_plan(args, run_data)

    # hmmpress the profile sets once for the whole run rather than re-parsing the HMM
    # text for every genome. The pressed files are scratch, so they live in a temp dir
    # that is cleaned up when the stage ends
    with tempfile.TemporaryDirectory(prefix="gtt-press-") as press_dir:
        if plan.do_scg:
            plan.pressed_scg_base = press_profiles(
                run_data.hmm_path, press_dir, "scg-profiles")
        if plan.do_pfam:
            plan.pressed_pfam_base = press_profiles(
                run_data.all_pfam_targets_hmm_path, press_dir, "pfam-profiles")

        run_data = _process_ncbi(args, run_data, plan)
        run_data = _process_genbank(args, run_data, plan)
        run_data = _process_fasta(args, run_data, plan)
        run_data = _process_amino_acid(args, run_data, plan)

    run_data = _finalize(args, run_data, plan)

    return run_data


def _process_ncbi(args, run_data, plan):
    if not run_data.ncbi_accs:
        return run_data

    from gtotree.utils.ncbi.parse_ncbi_assembly_summary import parse_assembly_summary
    from gtotree.utils.ncbi.get_ncbi_assembly_data import get_ncbi_assembly_summary_tab

    phase_stats.begin("processing genomes: ncbi accessions")
    report_processing_stage("ncbi", run_data)

    run_data = parse_assembly_summary(get_ncbi_assembly_summary_tab(), run_data)
    phase_stats.checkpoint("ncbi: after parsing assembly summary")

    to_process = [gd for gd in genomes_needing_processing(run_data.ncbi_accs, plan)
                  if gd.acc_was_found]

    if to_process:
        base_link_map = build_base_link_map(run_data)

        def preprocess(acc_gd, rd):
            return _process_one_ncbi_accession(acc_gd, rd, base_link_map)

        worker, apply_result = _fused(preprocess, _apply_ncbi_accession_status, plan)
        run_data = run_pooled_stage(to_process, worker, apply_result, args, run_data)
        run_data = capture_ncbi_failed_downloads(run_data)
        write_run_data(run_data)

    report_ncbi_update(run_data)
    return run_data


def _process_genbank(args, run_data, plan):
    if not args.genbank_files:
        return run_data

    phase_stats.begin("processing genomes: genbank files")
    report_processing_stage("genbank", run_data)
    to_process = genomes_needing_processing(run_data.genbank_files, plan)

    if to_process:
        worker, apply_result = _fused(_process_one_genbank_file, _apply_genbank_status, plan)
        run_data = run_pooled_stage(to_process, worker, apply_result, args, run_data)
        write_run_data(run_data)
        capture_failed_genbank_files(run_data)

    report_genbank_update(run_data)
    return run_data


def _process_fasta(args, run_data, plan):
    if not args.fasta_files:
        return run_data

    phase_stats.begin("processing genomes: fasta files")
    report_processing_stage("fasta", run_data)
    to_process = genomes_needing_processing(run_data.fasta_files, plan)

    if to_process:
        worker, apply_result = _fused(_process_one_fasta_file, _apply_fasta_status, plan)
        run_data = run_pooled_stage(to_process, worker, apply_result, args, run_data)
        write_run_data(run_data)
        capture_failed_fasta_files(run_data)

    report_fasta_update(run_data)
    return run_data


def _process_amino_acid(args, run_data, plan):
    if not args.amino_acid_files:
        return run_data

    phase_stats.begin("processing genomes: amino-acid files")
    report_processing_stage("amino-acid", run_data)
    to_process = genomes_needing_processing(run_data.amino_acid_files, plan)

    if to_process:
        worker, apply_result = _fused(_process_one_amino_acid_file,
                                      _apply_amino_acid_status, plan)
        run_data = run_pooled_stage(to_process, worker, apply_result, args, run_data)
        write_run_data(run_data)
        capture_failed_amino_acid_files(run_data)

    report_AA_update(run_data)
    return run_data


def _finalize(args, run_data, plan):
    """
    Build every combined output from the per-genome artifacts, then report
    """
    phase_stats.begin("processing genomes: combining outputs")
    report_processing_stage("processing-update", run_data)
    run_data.update_all_input_genomes()

    if plan.do_scg:
        capture_hmm_search_failures(run_data)
    report_genome_processing_update(run_data, searched=plan.do_scg)

    if plan.do_pfam:
        write_pfam_counts_table(run_data)
        print("") ; print("    Combining Pfam search results...")
        combine_all_pfam_hits(run_data.found_pfam_targets,
                              run_data.tmp_pfam_results_dir,
                              run_data.pfam_results_dir + "/pfam-hit-seqs")
    if run_data.target_pfams_file:
        run_data.additional_pfam_searching_done = True
        report_pfam_searching_update(run_data)

    if plan.do_ko:
        write_ko_counts_table(run_data)
        print("") ; print("    Combining KO search results...")
        combine_all_ko_hits(run_data.found_ko_targets,
                            run_data.tmp_ko_results_dir,
                            run_data.ko_results_dir + "/ko-hit-seqs")
    if run_data.target_kos_file:
        run_data.additional_ko_searching_done = True
        report_ko_searching_update(run_data)

    if plan.do_scg:
        phase_stats.checkpoint("combining: before SCG rebuild")
        run_data = rebuild_combined_SCG_outputs(run_data)
        phase_stats.checkpoint("combining: after SCG rebuild")

        run_data = check_target_SCGs_have_seqs(run_data, run_data.general_ext,
                                              SCGRemovalStage.NO_HITS)

    run_data.mark_stage_complete(PipelineStage.PROCESS_GENOMES)
    write_run_data(run_data)

    return run_data
