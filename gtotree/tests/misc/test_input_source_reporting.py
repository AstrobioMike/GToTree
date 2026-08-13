"""
Unit tests for the shared "where did these genomes come from" block in
gtotree/utils/misc/messaging.py.

The main GToTree driver's RUN INFO banner and `gtt search-pfams` / `gtt search-kos`
phase 1 both render this block. They used to describe the same inputs in two different
shapes ("- NCBI accessions listed in x.txt (8 genomes)" vs "8 NCBI accession(s) read
from x.txt"), which made the search subcommands read like a different program. Both now
go through these functions, so the two surfaces can't drift apart again.

The `-w` sub-bullets are read off RunData rather than off a RefGenomeSelection because
a resumed run doesn't re-resolve `-w` and has no selection object to ask.
"""

import argparse

from gtotree.utils.misc.general import (GenomeData, RunData, genome_source_label,
                                        genome_input_label,
                                        WANTED_REF_TAX_SOURCE_LABEL)
from gtotree.utils.misc.messaging import (input_genome_source_lines,
                                          total_input_genomes_line)


def _args(**overrides):
    base = {
        "wanted_ref_tax": None,
        "ncbi_accessions": None,
        "fasta_files": None,
        "amino_acid_files": None,
        "genbank_files": None,
    }
    base.update(overrides)
    return argparse.Namespace(**base)


class _FakeSelection:
    def __init__(self, accessions, canonical="Methanobacteria",
                 resolved_rank="class", effective_derep_rank="family"):
        self.accessions = accessions
        self.canonical = canonical
        self.resolved_rank = resolved_rank
        self.effective_derep_rank = effective_derep_rank
        self.warnings = []


def _run_data_with(accs=(), fastas=(), amino_acids=(), genbanks=()):
    rd = RunData()
    rd.ncbi_accs = [GenomeData.from_acc(a) for a in accs]
    rd.fasta_files = [GenomeData.from_path(p, "nt-fasta-file") for p in fastas]
    rd.amino_acid_files = [GenomeData.from_path(p, "aa-fasta-file")
                           for p in amino_acids]
    rd.genbank_files = [GenomeData.from_path(p, "genbank-file") for p in genbanks]
    rd.update_all_input_genomes()
    return rd


class TestInputGenomeSourceLines:

    def test_wanted_ref_tax_comes_first(self):
        """
        The `-w` block leads because it's the one input source the user didn't hand
        over as a file, so it's the one that needs explaining.
        """
        rd = _run_data_with(accs=["GCF_1", "GCF_2"])
        rd.merge_wanted_ref_tax_accessions(["GCA_9"])
        rd.record_wanted_ref_tax_selection(_FakeSelection(["GCA_9"]),
                                           taxon="Methanobacteria")

        lines = input_genome_source_lines(
            _args(wanted_ref_tax="Methanobacteria", ncbi_accessions="accs.txt"), rd)

        assert "--wanted-ref-tax 'Methanobacteria' (1 genomes)" in lines[0]
        assert "NCBI accessions listed in accs.txt (2 genomes)" in lines[2]

    def test_sources_keep_a_stable_order(self):
        rd = _run_data_with(accs=["GCF_1"], fastas=["a.fa"],
                            amino_acids=["b.faa"], genbanks=["c.gbk"])

        lines = input_genome_source_lines(
            _args(ncbi_accessions="accs.txt", fasta_files="f.txt",
                  amino_acid_files="aa.txt", genbank_files="gb.txt"), rd)

        assert [l.split(" listed in ")[0].strip() for l in lines] == [
            "- NCBI accessions",
            "- Fasta files",
            "- Amino-acid files",
            "- Genbank files",
        ]

    def test_derep_rank_is_reported_under_the_wanted_ref_tax_line(self):
        rd = _run_data_with()
        rd.merge_wanted_ref_tax_accessions(["GCA_1", "GCA_2"])
        rd.record_wanted_ref_tax_selection(
            _FakeSelection(["GCA_1", "GCA_2"], resolved_rank="class",
                           effective_derep_rank="family"))

        lines = input_genome_source_lines(_args(wanted_ref_tax="Methanobacteria"), rd)

        assert lines[1].strip() == ('- wanted rank "class" was dereplicated to one '
                                    'genome per family')

    def test_derep_off_says_so_rather_than_going_quiet(self):
        """
        Silence here is ambiguous: it reads the same as dereplication having happened
        and being unremarkable, when it's the case most likely to explain a surprising
        genome count.
        """
        rd = _run_data_with()
        rd.merge_wanted_ref_tax_accessions(["GCA_1"])
        rd.record_wanted_ref_tax_selection(
            _FakeSelection(["GCA_1"], resolved_rank="species",
                           effective_derep_rank=None))

        lines = input_genome_source_lines(_args(wanted_ref_tax="E. coli"), rd)

        assert "dereplication off" in lines[1]

    def test_accessions_already_listed_by_the_user_are_accounted_for(self):
        """
        A `-w` accession the user also listed with `-a` stays counted as theirs, so
        without this line the `-w` count silently disagrees with what the taxonomy
        selection reported.
        """
        rd = _run_data_with(accs=["GCA_1"])
        rd.merge_wanted_ref_tax_accessions(["GCA_1", "GCA_2"])
        rd.record_wanted_ref_tax_selection(_FakeSelection(["GCA_1", "GCA_2"]))

        lines = input_genome_source_lines(
            _args(wanted_ref_tax="Methanobacteria", ncbi_accessions="accs.txt"), rd)

        assert "(1 genomes)" in lines[0]
        assert "1 more were selected, but were already listed" in lines[2]

    def test_counts_add_up_to_the_reported_total(self):
        rd = _run_data_with(accs=["GCF_1", "GCF_2"], fastas=["a.fa"])
        rd.merge_wanted_ref_tax_accessions(["GCA_9"])
        rd.record_wanted_ref_tax_selection(_FakeSelection(["GCA_9"]))

        assert "Total input genomes: 4" in total_input_genomes_line(rd)

    def test_unused_sources_are_left_out(self):
        rd = _run_data_with(accs=["GCF_1"])

        lines = input_genome_source_lines(_args(ncbi_accessions="accs.txt"), rd)

        assert len(lines) == 1


class TestGenomeSourceLabels:

    def test_wanted_ref_tax_genomes_are_distinguishable_from_listed_accessions(self):
        """
        Both are accessions and both carry source == "accession", so a table keyed on
        the raw source alone can't answer "which of these did I actually ask for?".
        """
        rd = _run_data_with(accs=["GCF_1"])
        rd.merge_wanted_ref_tax_accessions(["GCA_9"])
        rd.record_wanted_ref_tax_selection(_FakeSelection(["GCA_9"]),
                                           taxon="Methanobacteria")

        listed, from_w = rd.ncbi_accs

        assert genome_source_label(listed) == "ncbi-accession"
        assert genome_source_label(from_w) == WANTED_REF_TAX_SOURCE_LABEL

    def test_wanted_ref_tax_genomes_report_the_taxon_as_their_input(self):
        rd = _run_data_with()
        rd.merge_wanted_ref_tax_accessions(["GCA_9"])
        rd.record_wanted_ref_tax_selection(_FakeSelection(["GCA_9"]),
                                           taxon="Methanobacteria")

        assert genome_input_label(rd.ncbi_accs[0], rd) == "-w Methanobacteria"

    def test_file_inputs_report_the_path_the_user_gave(self):
        rd = _run_data_with(fastas=["some/dir/genome_1.fa"])

        assert genome_input_label(rd.fasta_files[0], rd) == "some/dir/genome_1.fa"
