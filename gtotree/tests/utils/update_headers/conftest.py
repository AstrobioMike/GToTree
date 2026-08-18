"""
A completed GToTree run, on disk.

`gtt update-headers` reads a finished output directory rather than being handed a
RunData, so nearly every test here needs a real directory laid out the way a run
leaves it: run-files/run-data.json, the concatenated alignment in whichever of the two
places it belongs, a newick tree, and the summary table. Building that once here keeps
the tests about behavior rather than about fixture assembly.
"""

import os
import pytest  # type: ignore

from gtotree.utils.misc.general import GenomeData, RunData, write_run_data
from gtotree.utils.misc.stages import PipelineStage, GenomeRemovalStage


ACCESSIONS = ["GCF_000005845.2", "GCA_000009065.1", "GCF_000006945.2"]
DROPPED_ACCESSION = "GCF_000999999.1"

# a short alignment per genome; the actual residues are irrelevant here, only that
# each record's header round-trips
TAXIDS = {
    "GCF_000005845.2": "511145",
    "GCA_000009065.1": "83333",
    "GCF_000006945.2": "99287",
}

SEQS = {
    "GCF_000005845.2": "MKVLAAIL",
    "GCA_000009065.1": "MKVLAAVL",
    "GCF_000006945.2": "MKVLAAIM",
}


def build_completed_run(root, *, name="gtotree-output", mapping=None,
                        swapped=True, with_tree=True, alignment_in_run_files=None,
                        with_summary=True, nucleotide_mode=False,
                        target_pfams_file=None):
    """
    Write a finished run to `root/name` and return its path

    `mapping` is the label mapping the run itself produced ({} for a run that never
    re-labelled). `swapped` says whether that run's UPDATE_HEADERS stage completed --
    the two are separable, because a run whose swap failed still has a populated
    mapping dict but original IDs in its outputs.
    """
    out_dir = os.path.join(str(root), name)
    run_files = os.path.join(out_dir, "run-files")
    os.makedirs(run_files, exist_ok=True)

    ext = ".fasta" if nucleotide_mode else ".faa"

    run_data = RunData()
    run_data.ncbi_accs = []
    for acc in ACCESSIONS:
        gd = GenomeData.from_acc(acc)
        # set by parse_assembly_summary during the original run, and serialized into
        # run-data.json with the rest of the GenomeData
        gd.taxid = TAXIDS[acc]
        run_data.ncbi_accs.append(gd)

    dropped = GenomeData.from_acc(DROPPED_ACCESSION)
    dropped.mark_removed("test-dropped", GenomeRemovalStage.SCG_HIT_FILTER)
    run_data.ncbi_accs.append(dropped)

    run_data.update_all_input_genomes()

    run_data.output_dir = out_dir
    run_data.output_dir_rel = out_dir
    run_data.run_files_dir = run_files
    run_data.run_files_dir_rel = run_files
    run_data.run_data_path = os.path.join(run_files, "run-data.json")
    run_data.general_ext = ext
    run_data.nucleotide_mode = nucleotide_mode
    run_data.mapping_dict = dict(mapping or {})
    run_data.updating_headers = bool(mapping)
    run_data.target_pfams_file = target_pfams_file

    # a finished run's working dir is gone, and the NCBI sub-table with it; the
    # recorded path is left dangling exactly as a real finished run leaves it
    run_data.ncbi_sub_table_path = os.path.join(out_dir, "gtt-tmp-gone",
                                                "ncbi-accessions-info.tsv")

    for stage in (PipelineStage.PROCESS_GENOMES, PipelineStage.FILTER_GENES,
                  PipelineStage.FILTER_GENOMES, PipelineStage.ALIGN_SCG_SETS,
                  PipelineStage.CONCATENATE_SCG_SETS):
        run_data.mark_stage_complete(stage)

    labels_in_outputs = {acc: acc for acc in ACCESSIONS}
    if mapping and swapped:
        run_data.mark_stage_complete(PipelineStage.UPDATE_HEADERS)
        labels_in_outputs = {acc: mapping.get(acc, acc) for acc in ACCESSIONS}
    elif mapping:
        run_data.header_update_error = "something went wrong"

    # the original-ID alignment sits in the output dir until a swap moves it to
    # run-files/, and the relabeled copy takes its place
    if alignment_in_run_files is None:
        alignment_in_run_files = bool(mapping and swapped)

    original_dir = run_files if alignment_in_run_files else out_dir
    original_path = os.path.join(original_dir, f"aligned-SCGs{ext}")
    _write_fasta(original_path, {acc: SEQS[acc] for acc in ACCESSIONS})
    run_data.concatenated_alignment_path = original_path

    if alignment_in_run_files:
        final_alignment = os.path.join(out_dir, f"aligned-SCGs-mod-names{ext}")
        _write_fasta(final_alignment,
                     {labels_in_outputs[acc]: SEQS[acc] for acc in ACCESSIONS})
    else:
        final_alignment = original_path
    run_data.final_alignment_path = final_alignment

    if with_tree:
        tree_path = os.path.join(out_dir, f"{name}.tre")
        tips = [labels_in_outputs[acc] for acc in ACCESSIONS]
        with open(tree_path, "w") as f:
            f.write(f"(({tips[0]}:0.1,{tips[1]}:0.2)0.98:0.3,{tips[2]}:0.4);\n")
        run_data.final_tree_path = tree_path
        run_data.mark_stage_complete(PipelineStage.BUILD_TREE)

    run_data.mark_stage_complete(PipelineStage.FINALIZE)

    write_run_data(run_data)

    if with_summary:
        _write_summary(out_dir, labels_in_outputs)

    return out_dir


def _write_fasta(path, records):
    with open(path, "w") as f:
        for header, seq in records.items():
            f.write(f">{header}\n{seq}\n")


def _write_summary(out_dir, labels_in_outputs):
    """
    Enough of the run's genomes-summary-info.tsv to stand in for the real one
    """
    taxids = TAXIDS

    lines = ["genome_id\tinput\tsource\tlabel\ttaxid\tin_final_tree"]
    for acc in ACCESSIONS:
        lines.append("\t".join([acc, acc, "ncbi-accession",
                                labels_in_outputs[acc], taxids[acc], "Yes"]))
    lines.append("\t".join([DROPPED_ACCESSION, DROPPED_ACCESSION, "ncbi-accession",
                            DROPPED_ACCESSION, "NA", "No"]))

    with open(os.path.join(out_dir, "genomes-summary-info.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


@pytest.fixture
def completed_run(tmp_path):
    """A finished run that never re-labelled anything."""
    return build_completed_run(tmp_path)


@pytest.fixture
def relabeled_run(tmp_path):
    """A finished run that already swapped in labels from a mapping file."""
    return build_completed_run(tmp_path, mapping={
        "GCF_000005845.2": "E_coli_K12",
        "GCA_000009065.1": "E_coli_Sakai",
        "GCF_000006945.2": "S_enterica_LT2",
    })
