"""
The engine behind `gtt update-headers`.

A finished GToTree run already holds everything needed to re-label it: which genomes
went in, which of them made it into the tree, and what labels the outputs currently
carry. All of that is in `run-files/run-data.json`, so this program re-reads that
rather than re-deriving anything, and then does exactly what the main run's
`update_headers` stage does -- build a mapping dict from `-m` and/or a taxonomy
source, and swap it in.

Two things make this more than just calling that stage again:

  * the alignment can be relabeled from the original genome IDs, but the tree can't.
    There is no original-ID tree: whichever labels the alignment carried when the tree
    was built are the labels in the tree. So tree tips are re-labeled by going
    old-label -> genome id -> new label, using the mapping the previous run recorded.
  * every path in run-data.json is absolute, and a finished output directory is a
    thing people move and rename. `rehome_run_data` recomputes the directory fields
    from wherever the run is now, so a moved or renamed run still works.

Nothing here writes to the input directory. Library-raises / CLI-translates as
everywhere else: this module raises `UpdateHeadersError` and never prints or exits.
"""

import os
import glob

from gtotree.utils.misc.general import read_run_data, CorruptRunData, atomic_write_text
from gtotree.utils.misc.seqs import relabel_fasta, concatenated_alignment_locations
from gtotree.utils.misc.stages import PipelineStage
from gtotree.utils.misc.summary_info import generate_primary_summary_table


class UpdateHeadersError(Exception):
    """
    A user-facing problem, translated to a friendly message and a clean exit by the
    CLI layer.
    """


RUN_FILES_DIR_NAME = "run-files"
RUN_DATA_FILENAME = "run-data.json"
SUMMARY_FILENAME = "genomes-summary-info.tsv"
LABEL_MAP_FILENAME = "label-map.tsv"
ALIGNMENT_BASENAME = "aligned-SCGs-mod-names"

LABEL_MAP_COLUMNS = ("genome_id", "previous_label", "new_label")


################################################################################
# reading the completed run
################################################################################

def run_data_path_for(input_dir):
    return os.path.join(input_dir, RUN_FILES_DIR_NAME, RUN_DATA_FILENAME)


def load_completed_run(input_dir):
    """
    Read the RunData out of a finished GToTree output directory

    The three failure modes are kept distinct on purpose, because the fix is
    different for each: not a GToTree output directory at all, a run that never got
    far enough to have an alignment, and a run whose state file is damaged.
    """
    if not os.path.isdir(input_dir):
        raise UpdateHeadersError(
            f'The input directory "{input_dir}" (passed to `-i`) doesn\'t exist.')

    path = run_data_path_for(input_dir)
    if not os.path.isfile(path):
        raise UpdateHeadersError(
            f'"{input_dir}" (passed to `-i`) doesn\'t look like a GToTree output '
            f"directory -- there's no {RUN_FILES_DIR_NAME}/{RUN_DATA_FILENAME} in it. "
            "This wants the output directory of a completed GToTree run.")

    try:
        run_data = read_run_data(path)
    except CorruptRunData as e:
        raise UpdateHeadersError(
            f"We can't re-label that run because {e}. That usually means the run was "
            "interrupted while saving its state, in which case there's nothing here "
            "to work from.") from e

    if run_data is None:  # pragma: no cover - read_run_data only returns None if absent
        raise UpdateHeadersError(
            f'The run state at "{path}" couldn\'t be read.')

    if not run_data.stage_is_complete(PipelineStage.CONCATENATE_SCG_SETS):
        raise UpdateHeadersError(
            f'The GToTree run in "{input_dir}" didn\'t get as far as building its '
            "concatenated alignment, so there's nothing to re-label yet. Finish that "
            "run first (`GToTree ... -R`), then come back to this.")

    return run_data


def rehome_run_data(run_data, input_dir):
    """
    Point a loaded RunData's directory fields at where the run actually is now

    Every path in run-data.json is absolute and was recorded wherever the run was
    originally made, so a directory that's since been moved, renamed, or copied off
    a cluster would otherwise send us looking in a place that no longer exists.

    Only what this program reads gets re-homed. The working directory is deliberately
    blanked rather than re-pointed: a finished run deletes it, so anything that lived
    there is genuinely gone and pretending otherwise would just fail later and
    further away.
    """
    out_dir = os.path.abspath(input_dir.rstrip("/"))
    out_dir_rel = input_dir.rstrip("/")

    run_data.output_dir = out_dir
    run_data.output_dir_rel = out_dir_rel
    run_data.run_files_dir = os.path.join(out_dir, RUN_FILES_DIR_NAME)
    run_data.run_files_dir_rel = os.path.join(out_dir_rel, RUN_FILES_DIR_NAME)
    run_data.logs_dir = os.path.join(out_dir, "logs")
    run_data.individual_gene_alignments_dir = os.path.join(
        run_data.run_files_dir, "individual-alignments")
    run_data.pfam_results_dir = os.path.join(out_dir, "pfam-search-results")
    run_data.ko_results_dir = os.path.join(out_dir, "ko-search-results")

    run_data.tmp_dir = ""
    run_data.ncbi_sub_table_path = ""

    run_data.final_alignment_path = _relocate(run_data.final_alignment_path, out_dir)
    run_data.final_tree_path = _find_final_tree(run_data, out_dir)

    return run_data


def _is_usable(path):
    """
    Exists and is non-empty

    Deliberately not `file_is_usable_else_clear`, which deletes empty files: this
    program is read-only with respect to the run it's given.
    """
    try:
        return bool(path) and os.path.getsize(path) > 0
    except OSError:
        return False


def _relocate(recorded_path, *dirs):
    """
    The recorded file, found by basename in whichever of `dirs` actually holds it
    """
    if not recorded_path:
        return ""
    basename = os.path.basename(recorded_path)
    for directory in dirs:
        candidate = os.path.join(directory, basename)
        if _is_usable(candidate):
            return candidate
    return ""


def _find_final_tree(run_data, out_dir):
    """
    The finished tree, if there is one

    A run given `-N` has none at all, so an empty answer is a normal outcome rather
    than a problem. The tree is named after the output directory, so a *renamed*
    directory won't match by basename -- hence the fallback, which only commits when
    there's exactly one .tre to be unambiguous about.
    """
    found = _relocate(run_data.final_tree_path, out_dir)
    if found:
        return found

    candidates = [p for p in sorted(glob.glob(os.path.join(out_dir, "*.tre")))
                  if _is_usable(p)]
    if len(candidates) == 1:
        return candidates[0]
    return ""


def original_id_alignment_path(run_data):
    """
    The concatenated alignment as it was written with original genome IDs, or ""

    It sits in the output directory on a run that never swapped labels, and moves to
    run-files/ on one that did.
    """
    for path in concatenated_alignment_locations(run_data):
        if _is_usable(path):
            return path
    return ""


################################################################################
# labels
################################################################################

def labels_in_completed_outputs(run_data):
    """
    {genome id -> the label that genome carries in the finished tree and alignment}

    A run that never swapped labels, and one whose swap failed and fell back to the
    original IDs, both leave the outputs keyed by genome id -- and the second still
    has a populated `mapping_dict`, so the completion of the swap stage is what
    settles it, not the presence of a mapping.
    """
    swapped = (run_data.updating_headers
               and run_data.stage_is_complete(PipelineStage.UPDATE_HEADERS))

    if not swapped:
        return {gd.id: gd.id for gd in run_data.all_input_genomes}

    return {gd.id: run_data.mapping_dict.get(gd.id, gd.id)
            for gd in run_data.all_input_genomes}


def new_labels_for(run_data):
    """
    {genome id -> the label it should now carry}

    Genomes with no entry keep their genome id, which is what makes running this with
    neither `-m` nor a taxonomy flag a way to put the original IDs back.
    """
    return {gd.id: run_data.mapping_dict.get(gd.id, gd.id)
            for gd in run_data.all_input_genomes}


################################################################################
# writing the updated outputs
################################################################################

def output_name(prefix, basename):
    return f"{prefix}-{basename}" if prefix else basename


def write_updated_alignment(run_data, previous_labels, new_labels, out_path):
    """
    Write the re-labeled alignment, returning (source_path, num_records, num_relabeled)

    Preferred source is the original-ID concatenated alignment, because going
    straight from genome ids can't be thrown off by whatever the previous run
    labelled things. If run-files/ is gone, the finished alignment is relabeled
    instead, the same way the tree has to be.
    """
    source = original_id_alignment_path(run_data)

    if source:
        def label_for(header):
            return new_labels.get(header, header)
    else:
        source = run_data.final_alignment_path
        if not _is_usable(source):
            raise UpdateHeadersError(
                "We couldn't find the concatenated alignment from that run, in either "
                f'"{run_data.output_dir}" or "{run_data.run_files_dir}". If files were '
                "moved or removed since the run finished, there's nothing here to "
                "re-label.")
        label_for = _relabel_via_previous(previous_labels, new_labels)

    total, changed = relabel_fasta(source, out_path, label_for)

    return source, total, changed


def _relabel_via_previous(previous_labels, new_labels):
    """
    old label -> genome id -> new label

    Labels are unique within a run (both the mapping file and the taxonomy path
    guarantee it), so inverting the previous mapping is unambiguous.
    """
    id_by_previous = {label: gid for gid, label in previous_labels.items()}

    def label_for(header):
        genome_id = id_by_previous.get(header)
        if genome_id is None:
            return header
        return new_labels.get(genome_id, genome_id)

    return label_for


def write_updated_tree(run_data, previous_labels, new_labels, out_path):
    """
    Write the re-labeled tree, returning (num_tips, num_relabeled), or None if the
    run has no tree (`-N`)

    dendropy is used rather than text substitution so a label that happens to be a
    substring of another can't corrupt the newick, and it's read and written with the
    same underscore settings `gtt midpoint-root-tree` uses so the labels survive the
    round trip verbatim. Only leaf taxa are touched; support values live on nodes,
    not in the taxon namespace, so they pass through untouched.
    """
    import dendropy  # type: ignore

    source = run_data.final_tree_path
    if not _is_usable(source):
        return None

    tree = dendropy.Tree.get(path=source, schema="newick", preserve_underscores=True)

    label_for = _relabel_via_previous(previous_labels, new_labels)

    total = 0
    changed = 0
    for taxon in tree.taxon_namespace:
        if taxon.label is None:
            continue
        total += 1
        new_label = label_for(taxon.label)
        if new_label != taxon.label:
            taxon.label = new_label
            changed += 1

    tree.write(path=out_path, schema="newick", suppress_rooting=True,
               unquoted_underscores=True)

    return total, changed


def write_updated_summary_table(args, run_data, out_path):
    """
    The same genomes-summary-info.tsv a run produces, with the new labels in it

    The taxids ride along on the GenomeData loaded out of run-data.json, so this needs
    nothing from the original run's output files.
    """
    generate_primary_summary_table(args, run_data, out_path=out_path)


def write_label_map(run_data, previous_labels, new_labels, out_path):
    """
    A record of what each genome was called before and after

    Covers every input genome, not just the ones in the tree, so a genome that was
    named in the mapping file but dropped during the original run is visibly
    accounted for rather than silently missing.
    """
    def _write(f):
        f.write("\t".join(LABEL_MAP_COLUMNS) + "\n")
        for gd in run_data.all_input_genomes:
            f.write("\t".join([
                gd.id,
                previous_labels.get(gd.id, gd.id),
                new_labels.get(gd.id, gd.id),
            ]) + "\n")

    atomic_write_text(out_path, _write)


def regenerate_itol_files(run_data, summary_path, out_dir):
    """
    Re-make the iToL files for whichever additional searches the original run did,
    with the new labels, returning the directories written

    The iToL files are built from the hit-counts table plus the summary table, and
    it's the summary that carries the labels -- so pointing the same generator at the
    new summary and the original counts is all this takes.
    """
    from gtotree.utils.misc.itol import generate_search_itol_files

    written = []

    sources = []
    if run_data.target_pfams_file:
        sources.append((os.path.join(run_data.pfam_results_dir, "pfam-hit-counts.tsv"),
                        "pfam-search-results"))
    if run_data.target_kos_file:
        sources.append((os.path.join(run_data.ko_results_dir, "ko-hit-counts.tsv"),
                        "ko-search-results"))

    for counts_path, subdir in sources:
        if not _is_usable(counts_path):
            continue
        itol_dir = os.path.join(out_dir, subdir, "iToL-files")
        targets = generate_search_itol_files(counts_path=counts_path,
                                             summary_path=summary_path,
                                             itol_dir=itol_dir)
        if targets:
            written.append(itol_dir)

    return written
