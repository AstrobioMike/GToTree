"""
Output writing for `gtt gen-scg-hmms`

All writes go through `_atomic_write`, so an interrupted run never leaves a truncated
table

Files produced in the output directory:

    <name>.hmm                      the new SCG-HMM set (the actual deliverable)
    SCG-targets-info.tsv            one row per retained target: acc, name, description,
                                    coverage, and its single-/multi-copy genome counts
    Pfam-hit-counts.tsv             genome x Pfam hit-count matrix
    target-genomes.tsv              the genomes actually used, with their source
    shared-protein-pairs.tsv        every pair of searched Pfams that hit one protein
                                    together, and in how many genomes
    shared-protein-exclusions.tsv   Pfams that passed the single-copy cutoff but were
                                    dropped for commonly hitting the same protein as a
                                    retained one (header only if there were none)
    selection-params.tsv            the parameters the targets were selected with
    pfam-version-used.txt           the Pfam release the set was built from
    removed-genomes.tsv             any input genome that dropped out, and where/why

That last one isn't written here: it comes from the shared
`summary_info.write_removed_genomes_report`, so this program's account of what it lost
has the same name, columns, and wording as the main driver's and the search
subcommands'.

The scan-record tables (`SCAN_RECORD_FILENAMES`) are everything `--from-run` needs to
select a new set without redoing the search, and are copied as-is into its output dir,
so a `--from-run` output dir is itself a complete run that can be selected from again.
`SCG-targets-info.tsv` is written last in both cases: its presence is what marks a run
as finished.
"""

import os
import shutil
import contextlib

from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import (
    GenSCGHMMsError,
    shared_pair_key,
    split_shared_pair_key,
)
from gtotree.utils.misc.messaging import REMOVED_GENOMES_FILENAME


HMM_INFO_FILENAME = "SCG-targets-info.tsv"
HIT_COUNTS_FILENAME = "Pfam-hit-counts.tsv"
TARGET_GENOMES_FILENAME = "target-genomes.tsv"
PFAM_VERSION_FILENAME = "pfam-version-used.txt"
SHARED_PROTEIN_EXCLUSIONS_FILENAME = "shared-protein-exclusions.tsv"
SHARED_PROTEIN_PAIRS_FILENAME = "shared-protein-pairs.tsv"
SELECTION_PARAMS_FILENAME = "selection-params.tsv"

# what a finished run needs to have for `--from-run` to select from it
SCAN_RECORD_FILENAMES = (
    HIT_COUNTS_FILENAME,
    TARGET_GENOMES_FILENAME,
    SHARED_PROTEIN_PAIRS_FILENAME,
    PFAM_VERSION_FILENAME,
)


@contextlib.contextmanager
def _atomic_write(path):
    """ Write to `path` via a `.part` temp, moving into place only on success. """
    tmp_path = path + ".part"
    handle = open(tmp_path, "w")
    try:
        yield handle
        handle.close()
        os.replace(tmp_path, path)
    except BaseException:
        handle.close()
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise


def write_scg_targets_info(out_dir, wanted_accs, pfam_info, single_copy_counts,
                           multi_copy_counts, num_genomes):
    """
    Write the info table for the Pfams retained as SCG targets, with the number (and
    percent) of genomes each was single-copy in, which is what they were selected and
    ranked on, plus the number of genomes each was multi-copy in.
    """
    path = os.path.join(out_dir, HMM_INFO_FILENAME)
    with _atomic_write(path) as out:
        out.write("pfam_id\tname\tdescription\taverage_coverage\t"
                  "num_genomes_single_copy\tperc_genomes_single_copy\t"
                  "num_genomes_multi_copy\n")
        for acc in wanted_accs:
            info = pfam_info.get(acc)
            if info is None:
                name, description, coverage = "NA", "NA", "NA"
            else:
                name, description, coverage = info.name, info.description, info.coverage
            n_single = single_copy_counts.get(acc, 0)
            perc = round(n_single / num_genomes * 100, 2) if num_genomes else 0
            out.write("\t".join([
                acc, name, description, str(coverage),
                str(n_single), f"{perc:g}", str(multi_copy_counts.get(acc, 0)),
            ]) + "\n")
    return path


def write_hit_counts(out_dir, genome_ids, filtered_accs, per_genome_counts):
    """
    Write the full genome x Pfam hit-count matrix.

    Rows are genomes, columns are every profile that was searched (not just the
    retained ones), so the table can be used to see why a marker was dropped.
    """
    path = os.path.join(out_dir, HIT_COUNTS_FILENAME)
    with _atomic_write(path) as out:
        out.write("genome\t" + "\t".join(filtered_accs) + "\n")
        for genome_id in genome_ids:
            counts = per_genome_counts.get(genome_id, {})
            row = "\t".join(str(counts.get(acc, 0)) for acc in filtered_accs)
            out.write(f"{genome_id}\t{row}\n")
    return path


def read_hit_counts(path):
    """
    Read back a `Pfam-hit-counts.tsv`.

    Returns (genome_ids, searched_accs, hits_by_genome), with hits_by_genome in the
    same {genome_id: {acc: count}} shape the search produces (zero counts left out),
    so it can go straight into `select_scg_targets`.
    """
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if len(header) < 2 or header[0] != "genome":
            raise GenSCGHMMsError(
                f"'{path}' doesn't look like a {HIT_COUNTS_FILENAME} table.")
        accs = header[1:]

        genome_ids, hits_by_genome = [], {}
        for line_num, line in enumerate(f, 2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != len(header):
                raise GenSCGHMMsError(
                    f"line {line_num} of '{path}' has {len(parts)} columns, but the "
                    f"header has {len(header)}.")
            try:
                counts = [int(x) for x in parts[1:]]
            except ValueError:
                raise GenSCGHMMsError(
                    f"line {line_num} of '{path}' holds a non-integer count.") from None
            genome_ids.append(parts[0])
            hits_by_genome[parts[0]] = {acc: n for acc, n in zip(accs, counts) if n}

    if not genome_ids:
        raise GenSCGHMMsError(f"'{path}' holds no genomes.")

    return genome_ids, accs, hits_by_genome


def write_shared_protein_pairs(out_dir, genomes_sharing, num_genomes, pfam_info):
    """
    Write every pair of searched Pfams that hit the same protein in at least one
    genome, and in how many. Most-shared first. Written even when empty.

    Unlike `shared-protein-exclusions.tsv`, this isn't limited to the Pfams that were
    retained, so it holds what's needed to resolve conflicts under any other cutoffs.
    """
    def _name(acc):
        info = pfam_info.get(acc)
        return info.name if info is not None else "NA"

    rows = sorted(((key, n) for key, n in genomes_sharing.items() if n > 0),
                  key=lambda item: (-item[1], item[0]))

    path = os.path.join(out_dir, SHARED_PROTEIN_PAIRS_FILENAME)
    with _atomic_write(path) as out:
        out.write("pfam_id_1\tname_1\tpfam_id_2\tname_2\t"
                  "num_genomes_with_both_in_one_protein\tperc_genomes_with_both_in_one_protein\n")
        for key, n in rows:
            a, b = split_shared_pair_key(key)
            perc = round(n / num_genomes * 100, 2) if num_genomes else 0
            out.write("\t".join([a, _name(a), b, _name(b), str(n), f"{perc:g}"]) + "\n")
    return path


def read_shared_protein_pairs(path):
    """
    Read back a `shared-protein-pairs.tsv` as {pair_key: num_genomes}, the shape
    `resolve_shared_protein_conflicts` takes.
    """
    genomes_sharing = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            a_col = header.index("pfam_id_1")
            b_col = header.index("pfam_id_2")
            n_col = header.index("num_genomes_with_both_in_one_protein")
        except ValueError:
            raise GenSCGHMMsError(
                f"'{path}' doesn't look like a {SHARED_PROTEIN_PAIRS_FILENAME} "
                "table.") from None
        for line_num, line in enumerate(f, 2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            try:
                n = int(parts[n_col])
                key = shared_pair_key(parts[a_col], parts[b_col])
            except (IndexError, ValueError):
                raise GenSCGHMMsError(
                    f"line {line_num} of '{path}' couldn't be read.") from None
            genomes_sharing[key] = n
    return genomes_sharing


def write_selection_params(out_dir, params):
    """
    Record the parameters the targets were selected with, as `parameter<TAB>value`
    rows from the (name, value) pairs in `params`. Once one scan can produce several
    sets via `--from-run`, this is what tells them apart.
    """
    path = os.path.join(out_dir, SELECTION_PARAMS_FILENAME)
    with _atomic_write(path) as out:
        out.write("parameter\tvalue\n")
        for name, value in params:
            out.write(f"{name}\t{'NA' if value is None else value}\n")
    return path


def check_previous_run(run_dir):
    """
    Make sure `run_dir` is a finished `gen-scg-hmms` run holding everything
    `--from-run` needs, raising GenSCGHMMsError naming what's missing if not.
    """
    if not os.path.isdir(run_dir):
        raise GenSCGHMMsError(f"The `--from-run` directory '{run_dir}' can't be found.")

    if not os.path.isfile(os.path.join(run_dir, HMM_INFO_FILENAME)):
        raise GenSCGHMMsError(
            f"'{run_dir}' doesn't hold a finished `gen-scg-hmms` run (there's no "
            f"{HMM_INFO_FILENAME}). If it was interrupted, it can be finished with "
            "`-R`/`--resume` first.")

    missing = [name for name in SCAN_RECORD_FILENAMES
               if not os.path.isfile(os.path.join(run_dir, name))]
    if missing:
        raise GenSCGHMMsError(
            f"'{run_dir}' is missing {', '.join(missing)}, which `--from-run` needs. "
            "It was likely made by an earlier version of `gen-scg-hmms`, so it would "
            "need to be re-run to be selected from.")


def copy_scan_record(run_dir, out_dir):
    """
    Copy the previous run's scan-record tables (and its removed-genomes report, if it
    has one) into `out_dir`, so the new output dir is a complete run of its own.
    """
    names = list(SCAN_RECORD_FILENAMES)
    if os.path.isfile(os.path.join(run_dir, REMOVED_GENOMES_FILENAME)):
        names.append(REMOVED_GENOMES_FILENAME)

    for name in names:
        dest = os.path.join(out_dir, name)
        tmp_path = dest + ".part"
        try:
            shutil.copyfile(os.path.join(run_dir, name), tmp_path)
            os.replace(tmp_path, dest)
        except BaseException:
            try:
                os.remove(tmp_path)
            except OSError:
                pass
            raise
    return names


def write_shared_protein_exclusions(out_dir, exclusions, pfam_info):
    """
    Write the Pfams dropped for sharing proteins with a retained Pfam, and which one
    they lost out to. The count columns are the number (and percent) of target genomes
    in which both Pfams hit the same single protein. Written even when empty.
    """
    def _name(acc):
        info = pfam_info.get(acc)
        return info.name if info is not None else "NA"

    path = os.path.join(out_dir, SHARED_PROTEIN_EXCLUSIONS_FILENAME)
    with _atomic_write(path) as out:
        out.write("dropped_pfam_id\tdropped_name\tkept_pfam_id\tkept_name\t"
                  "num_genomes_with_both_in_one_protein\tperc_genomes_with_both_in_one_protein\n")
        for ex in exclusions:
            out.write("\t".join([
                ex.dropped_acc, _name(ex.dropped_acc),
                ex.kept_acc, _name(ex.kept_acc),
                str(ex.num_genomes_shared), f"{ex.percent_genomes_shared:g}",
            ]) + "\n")
    return path


def write_target_genomes(out_dir, genome_ids, run_data):
    """
    Write the genomes that made it through and were actually searched
    """
    from gtotree.utils.misc.general import genome_input_label, genome_source_label

    path = os.path.join(out_dir, TARGET_GENOMES_FILENAME)
    by_id = {gd.id: gd for gd in run_data.all_input_genomes}

    with _atomic_write(path) as out:
        out.write("genome_id\tinput\tsource\torganism_name\n")
        for genome_id in genome_ids:
            gd = by_id.get(genome_id)
            if gd is None:
                out.write(f"{genome_id}\tNA\tNA\tNA\n")
                continue
            out.write("\t".join([
                gd.id,
                genome_input_label(gd, run_data),
                genome_source_label(gd),
                gd.organism_name or "NA",
            ]) + "\n")
    return path


def write_pfam_version(out_dir, version):
    """ Record the Pfam release the SCG set was built from. """
    path = os.path.join(out_dir, PFAM_VERSION_FILENAME)
    with _atomic_write(path) as out:
        out.write(f"{version}\n")
    return path


def default_hmm_filename(output_dir, num_targets):
    """
    Pick the output HMM filename
    """
    base = os.path.basename(os.path.normpath(output_dir))
    if base and base != "gtt-gen-scg-hmms-output":
        return f"{base}.hmm"
    return f"wanted-{num_targets}-scg-targets.hmm"
