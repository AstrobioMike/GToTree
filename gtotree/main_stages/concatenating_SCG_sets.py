import os

from gtotree.utils.misc.general import file_is_usable_else_clear, write_run_data
from gtotree.utils.misc.messaging import report_processing_stage
from gtotree.utils.misc.stages import PipelineStage
from gtotree.utils.misc.seqs import (concatenate_alignments,
                                gen_partitions_file,
                                get_alignment_length)


def concatenate_SCG_sets(run_data):

    report_processing_stage("concatenate-SCG-sets", run_data)

    if not _concatenation_can_be_skipped(run_data):

        run_data, dict_of_genomes, SCG_IDs = concatenate_alignments(run_data)

        run_data.final_alignment_length = get_alignment_length(run_data.concatenated_alignment_path)

        gen_partitions_file(run_data, SCG_IDs, dict_of_genomes)

        run_data.mark_stage_complete(PipelineStage.CONCATENATE_SCG_SETS)
        write_run_data(run_data)

    return run_data


def _concatenation_can_be_skipped(run_data):
    """
    True when a previous run finished concatenating and its output is still there
    """
    if not run_data.stage_is_complete(PipelineStage.CONCATENATE_SCG_SETS):
        return False
    return _concatenated_alignment_is_present(run_data)


def _concatenated_alignment_is_present(run_data):
    """
    The concatenated alignment legitimately lives at either of two paths
    """
    path = run_data.concatenated_alignment_path
    if not path:
        return False
    if file_is_usable_else_clear(path):
        return True
    if not run_data.run_files_dir:
        return False
    return file_is_usable_else_clear(
        os.path.join(run_data.run_files_dir, os.path.basename(path)))
