"""
Unit tests for the input genome-ID collision check in
gtotree/utils/misc/preflight_checks.py.

The genome ID is derived from the input filename and then used as both a filename and a
sequence-header prefix, so two inputs resolving to one ID silently clobber each other
during preprocessing while every count downstream still reports both. `check_for_duplicates`
doesn't cover this -- it only collapses identical lines within a single input file.
"""

import pytest  # type: ignore

from gtotree.utils.misc.general import GenomeData, RunData
from gtotree.utils.misc.preflight_checks import (check_for_genome_id_collisions,
                                                 collect_genome_id_collisions)


def _run_data(**by_source):
    """
    by_source: {source_field: [paths or accessions]}
    """
    rd = RunData()
    for field, entries in by_source.items():
        source = {"ncbi_accs": "ncbi-accession",
                  "genbank_files": "genbank-file",
                  "fasta_files": "nucleotide-fasta",
                  "amino_acid_files": "amino-acid-fasta"}[field]
        for entry in entries:
            gd = (GenomeData.from_acc(entry) if field == "ncbi_accs"
                  else GenomeData.from_path(entry, source))
            getattr(rd, field).append(gd)
    rd.update_all_input_genomes()
    return rd


def _collision_ids(rd):
    return sorted(collect_genome_id_collisions(rd))


class TestCollectGenomeIDCollisions:

    def test_distinct_names_do_not_collide(self):
        rd = _run_data(fasta_files=["a/genome_1.fa", "b/genome_2.fa"])
        assert collect_genome_id_collisions(rd) == {}

    def test_same_basename_in_different_directories(self):
        """
        The case `check_for_duplicates` can't see: the two lines differ, so it passes
        them both through, but they reduce to one ID.
        """
        rd = _run_data(fasta_files=["set_A/genome_1.fa", "set_B/genome_1.fa"])
        assert _collision_ids(rd) == ["genome_1"]

    def test_same_stem_with_different_sequence_extensions(self):
        rd = _run_data(fasta_files=["d/genome_1.fa", "d/genome_1.fasta"])
        assert _collision_ids(rd) == ["genome_1"]

    def test_gzipped_and_plain_versions_of_one_name(self):
        rd = _run_data(fasta_files=["d/genome_1.fa", "d/genome_1.fa.gz"])
        assert _collision_ids(rd) == ["genome_1"]

    def test_across_two_input_sources(self):
        """
        Per-file validation is structurally unable to catch this one, since the two
        entries arrive through different flags.
        """
        rd = _run_data(fasta_files=["d/genome_1.fa"],
                       amino_acid_files=["d/genome_1.faa"])
        assert _collision_ids(rd) == ["genome_1"]

    def test_a_file_named_after_an_accession_collides_with_it(self):
        rd = _run_data(ncbi_accs=["GCF_000005845.2"],
                       fasta_files=["d/GCF_000005845.2.fa"])
        assert _collision_ids(rd) == ["GCF_000005845.2"]

    def test_every_colliding_entry_is_collected_not_just_the_extras(self):
        rd = _run_data(fasta_files=["a/g.fa", "b/g.fa", "c/g.fa"])
        assert len(collect_genome_id_collisions(rd)["g"]) == 3

    def test_several_separate_collisions_are_reported_together(self):
        rd = _run_data(fasta_files=["a/g1.fa", "b/g1.fa", "a/g2.fa", "b/g2.fa"])
        assert _collision_ids(rd) == ["g1", "g2"]

    def test_a_clean_set_of_genomes_is_untouched(self):
        rd = _run_data(ncbi_accs=["GCF_000005845.2", "GCA_000008865.2"],
                       genbank_files=["x.gb"],
                       fasta_files=["y.fa"],
                       amino_acid_files=["z.faa"])
        assert collect_genome_id_collisions(rd) == {}


class TestCheckForGenomeIDCollisions:

    def test_a_clean_set_passes_silently(self, capsys):
        check_for_genome_id_collisions(_run_data(fasta_files=["a.fa", "b.fa"]))
        assert capsys.readouterr().out == ""

    def test_a_collision_stops_the_run(self):
        rd = _run_data(fasta_files=["set_A/genome_1.fa", "set_B/genome_1.fa"])
        with pytest.raises(SystemExit):
            check_for_genome_id_collisions(rd)

    def test_the_offending_paths_are_named(self, capsys):
        """
        The ID alone doesn't tell you which files to go rename -- that's the whole
        problem being reported -- so both provided paths have to appear.
        """
        rd = _run_data(fasta_files=["set_A/genome_1.fa", "set_B/genome_1.fa"])
        with pytest.raises(SystemExit):
            check_for_genome_id_collisions(rd)

        out = capsys.readouterr().out
        assert "set_A/genome_1.fa" in out
        assert "set_B/genome_1.fa" in out

    def test_each_offender_is_tagged_with_the_flag_it_came_in_on(self, capsys):
        rd = _run_data(fasta_files=["d/genome_1.fa"],
                       amino_acid_files=["d/genome_1.faa"])
        with pytest.raises(SystemExit):
            check_for_genome_id_collisions(rd)

        out = capsys.readouterr().out
        assert "-f  d/genome_1.fa" in out
        assert "-A  d/genome_1.faa" in out
