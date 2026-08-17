"""
Unit tests for the shared "where did these genomes come from" block in
gtotree/utils/misc/messaging.py.

The main GToTree driver's RUN INFO banner and `gtt search-pfams` / `gtt search-kos`
phase 1 both render this block. They used to describe the same inputs in two different
shapes ("- NCBI accessions listed in x.txt (8 genomes)" vs "8 NCBI accession(s) read
from x.txt"), which made the search subcommands read like a different program. Both now
go through these functions, so the two surfaces can't drift apart again.

The `-w` sub-bullets are read off RunData rather than off a RefGenomeSelection because
a resumed run doesn't re-resolve `-w` and has no selection object to ask. They are a
LIST of selections because `-w` is repeatable in gen-scg-hmms / search-pfams /
search-kos: each taxon is resolved and dereplicated on its own and gets its own line.
"""

import argparse
import os

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


def _add_wanted(rd, accessions, taxon="Methanobacteria", **selection_kwargs):
    """Merge + record one `-w` selection the way the CLIs do."""
    selection = _FakeSelection(accessions, canonical=taxon, **selection_kwargs)
    added = rd.merge_wanted_ref_tax_accessions(accessions, taxon=taxon)
    rd.record_wanted_ref_tax_selection(selection, taxon=taxon, num_added=added)
    return selection


def _run_data_with(accs=(), fastas=(), amino_acids=(), genbanks=()):
    rd = RunData()
    rd.ncbi_accs = [GenomeData.from_acc(a) for a in accs]
    rd.fasta_files = [GenomeData.from_path(p, "nucleotide-fasta") for p in fastas]
    rd.amino_acid_files = [GenomeData.from_path(p, "amino-acid-fasta")
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
        _add_wanted(rd, ["GCA_9"])

        lines = input_genome_source_lines(
            _args(wanted_ref_tax=["Methanobacteria"], ncbi_accessions="accs.txt"), rd)

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
        _add_wanted(rd, ["GCA_1", "GCA_2"], resolved_rank="class",
                    effective_derep_rank="family")

        lines = input_genome_source_lines(
            _args(wanted_ref_tax=["Methanobacteria"]), rd)

        assert lines[1].strip() == ('- input rank, class, was dereplicated to one '
                                    'genome per family')

    def test_derep_off_says_so_rather_than_going_quiet(self):
        """
        Silence here is ambiguous: it reads the same as dereplication having happened
        and being unremarkable, when it's the case most likely to explain a surprising
        genome count.
        """
        rd = _run_data_with()
        _add_wanted(rd, ["GCA_1"], taxon="E. coli", resolved_rank="species",
                    effective_derep_rank=None)

        lines = input_genome_source_lines(_args(wanted_ref_tax=["E. coli"]), rd)

        assert "--derep-rank off" in lines[1]

    def test_accessions_already_listed_by_the_user_are_accounted_for(self):
        """
        A `-w` accession the user also listed with `-a` stays counted as theirs, so
        without this line the `-w` count silently disagrees with what the taxonomy
        selection reported.
        """
        rd = _run_data_with(accs=["GCA_1"])
        _add_wanted(rd, ["GCA_1", "GCA_2"])

        lines = input_genome_source_lines(
            _args(wanted_ref_tax=["Methanobacteria"], ncbi_accessions="accs.txt"), rd)

        assert "(1 genomes)" in lines[0]
        assert "1 more were selected, but were already counted" in lines[2]

    def test_counts_add_up_to_the_reported_total(self):
        rd = _run_data_with(accs=["GCF_1", "GCF_2"], fastas=["a.fa"])
        _add_wanted(rd, ["GCA_9"])

        assert "Total input genomes: 4" in total_input_genomes_line(rd)

    def test_unused_sources_are_left_out(self):
        rd = _run_data_with(accs=["GCF_1"])

        lines = input_genome_source_lines(_args(ncbi_accessions="accs.txt"), rd)

        assert len(lines) == 1


class TestGenomeSourceLabels:

    def test_wanted_ref_tax_genomes_are_distinguishable_from_listed_accessions(self):
        """
        Both are accessions and both carry source == "ncbi-accession", so a table keyed on
        the raw source alone can't answer "which of these did I actually ask for?".
        """
        rd = _run_data_with(accs=["GCF_1"])
        _add_wanted(rd, ["GCA_9"])

        listed, from_w = rd.ncbi_accs

        assert genome_source_label(listed) == "ncbi-accession"
        assert genome_source_label(from_w) == WANTED_REF_TAX_SOURCE_LABEL

    def test_wanted_ref_tax_genomes_report_the_taxon_as_their_input(self):
        rd = _run_data_with()
        _add_wanted(rd, ["GCA_9"])

        assert genome_input_label(rd.ncbi_accs[0], rd) == "-w Methanobacteria"

    def test_each_wanted_ref_tax_taxon_reports_its_own_genomes(self):
        """
        `-w` is repeatable, so one run-level taxon can't say which request brought a
        given accession in. Each selection gets its own line, and each genome carries
        the taxon that selected it.
        """
        rd = _run_data_with()
        _add_wanted(rd, ["GCA_1", "GCA_2"], taxon="Bacteria")
        _add_wanted(rd, ["GCA_2", "GCA_3"], taxon="Archaea")

        lines = input_genome_source_lines(
            _args(wanted_ref_tax=["Bacteria", "Archaea"]), rd)

        assert "'Bacteria' (2 genomes)" in lines[0]
        assert "'Archaea' (1 genomes)" in lines[2]
        # GCA_2 was selected by both; the first taxon to select it owns it
        assert "1 more were selected, but were already counted" in lines[4]

        by_id = {gd.id: gd for gd in rd.ncbi_accs}
        assert genome_input_label(by_id["GCA_2"], rd) == "-w Bacteria"
        assert genome_input_label(by_id["GCA_3"], rd) == "-w Archaea"

    def test_file_inputs_report_the_path_the_user_gave(self):
        rd = _run_data_with(fastas=["some/dir/genome_1.fa"])

        assert genome_input_label(rd.fasta_files[0], rd) == "some/dir/genome_1.fa"


class TestPhaseOneIsSharedAcrossPrograms:
    """
    `gtt gen-scg-hmms`, `gtt search-pfams`, and `gtt search-kos` all describe their
    input genomes in phase 1. They used to do it in two different shapes, which made
    gen-scg-hmms read like a different program. They now share one implementation, and
    these check that they still do -- a second hand-rolled rendering is exactly the
    thing that drifts.
    """

    def _inputs(self, tmp_path):
        genome = tmp_path / "g1.faa"
        genome.write_text(">p1\nMKVLAAAL\n")
        listing = tmp_path / "aa.txt"
        listing.write_text(f"{genome}\n")
        return str(listing)

    def _args(self, listing):
        return argparse.Namespace(
            ncbi_accessions=None, wanted_ref_tax=None, genbank_files=None,
            fasta_files=None, amino_acid_files=listing, source="gtdb",
            target_rank=None, derep_rank="off",
            min_completeness=None, max_contamination=None)

    def test_both_programs_render_the_same_phase_one(self, tmp_path, capsys):
        from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_genomes import build_run_data
        from gtotree.utils.hmms.gen_scg_hmms import gen_scg_hmms_cli as gen_cli
        from gtotree.utils.target_search import target_search_stages as stages
        from gtotree.utils.target_search.target_search_setup import (
            _populate_input_genomes)

        listing = self._inputs(tmp_path)
        args = self._args(listing)

        out_dir = str(tmp_path / "out")
        work_dir = str(tmp_path / "work")
        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(work_dir, exist_ok=True)

        gen_cli.phase_resolve_genomes(args, build_run_data(args, out_dir, work_dir))
        gen_scg_hmms_output = capsys.readouterr().out

        search_run_data = RunData()
        _populate_input_genomes(args, search_run_data)
        stages.resolve_input_genomes(args, search_run_data)
        search_output = capsys.readouterr().out

        assert gen_scg_hmms_output == search_output

    def test_both_programs_use_the_same_resolver(self):
        """
        The identical-output test above would still pass if one program grew its own
        copy that happened to match today. This pins the actual sharing.
        """
        from gtotree.utils.misc.general import resolve_input_genomes
        from gtotree.utils.target_search import target_search_stages as stages
        import inspect

        assert "resolve_input_genomes" in inspect.getsource(stages.resolve_input_genomes)
        assert resolve_input_genomes.__module__ == "gtotree.utils.misc.general"
