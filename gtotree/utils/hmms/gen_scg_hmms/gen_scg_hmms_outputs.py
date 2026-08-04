"""
Output writing for `gtt gen-scg-hmms`

All writes go through `_atomic_write`, so an interrupted run never leaves a truncated
table

Files produced in the output directory:

    <name>.hmm                      the new SCG-HMM set (the actual deliverable)
    SCG-targets-info.tsv            one row per retained target: acc, name, description
    Pfam-hit-counts.tsv             genome x Pfam hit-count matrix
    missed-accessions.tsv           any input accession that dropped out, and why
    target-genomes.tsv              the genomes actually used, with their source
    pfam-version-used.txt           the Pfam release the set was built from
"""

import os
import contextlib


HMM_INFO_FILENAME = "SCG-targets-info.tsv"
HIT_COUNTS_FILENAME = "Pfam-hit-counts.tsv"
MISSED_FILENAME = "missed-accessions.tsv"
TARGET_GENOMES_FILENAME = "target-genomes.tsv"
PFAM_VERSION_FILENAME = "pfam-version-used.txt"


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


def write_scg_targets_info(out_dir, wanted_accs, pfam_info):
    """ Write the info table for the Pfams retained as SCG targets. """
    path = os.path.join(out_dir, HMM_INFO_FILENAME)
    with _atomic_write(path) as out:
        out.write("pfam_id\tname\tdescription\taverage_coverage\n")
        for acc in wanted_accs:
            info = pfam_info.get(acc)
            if info is None:
                out.write(f"{acc}\tNA\tNA\tNA\n")
            else:
                out.write(f"{acc}\t{info.name}\t{info.description}\t{info.coverage}\n")
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


def write_missed_accessions(out_dir, missed):
    """
    Write the accessions that didn't make it, with the reason for each.

    `missed` is an iterable of (accession, reason). Returns the path, or None if
    nothing was missed (in which case no file is written).
    """
    missed = list(missed)
    if not missed:
        return None

    path = os.path.join(out_dir, MISSED_FILENAME)
    with _atomic_write(path) as out:
        out.write("accession\twhy_missing\n")
        for acc, reason in missed:
            out.write(f"{acc}\t{reason}\n")
    return path


def write_target_genomes(out_dir, genome_ids, sources=None, organism_names=None):
    """ Write the genomes that made it through and were actually searched. """
    path = os.path.join(out_dir, TARGET_GENOMES_FILENAME)
    sources = sources or {}
    organism_names = organism_names or {}
    with _atomic_write(path) as out:
        out.write("accession\tsource\torganism_name\n")
        for genome_id in genome_ids:
            source = sources.get(genome_id, "NA")
            organism = organism_names.get(genome_id) or "NA"
            out.write(f"{genome_id}\t{source}\t{organism}\n")
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
