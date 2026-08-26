"""
`--exclusion-list` drops accessions from the references `-w` pulls in, before they
reach processing.

The exclusion is applied inside merge_wanted_ref_tax(), against each `-w` selection's
accessions only, and never against accessions the user provided directly through `-a`.
The file itself is validated in preflight exactly like an `-a` accessions file, and is
only meaningful alongside `-w`.
"""

import pytest  # type: ignore

from gtotree.utils.misc.general import RunData, GenomeData
from gtotree.utils.misc import preflight_checks
from gtotree.utils.misc.preflight_checks import (merge_wanted_ref_tax,
                                            check_input_genome_files,
                                            build_fingerprint)
from gtotree.utils.misc.messaging import _wanted_ref_tax_detail_lines


class _Selection:
    """Stand-in for a RefGenomeSelection as merge/record consume it."""
    def __init__(self, accessions, canonical="Alteromonas",
                 resolved_rank="genus", effective_derep_rank="species", warnings=()):
        self.accessions = list(accessions)
        self.canonical = canonical
        self.resolved_rank = resolved_rank
        self.effective_derep_rank = effective_derep_rank
        self.warnings = list(warnings)


class _Args:
    def __init__(self, **kw):
        self.wanted_ref_tax = None
        self.ncbi_accessions = None
        self.genbank_files = None
        self.fasta_files = None
        self.amino_acid_files = None
        self.exclusion_list = None
        self.__dict__.update(kw)


def _write(path, entries):
    path.write_text("\n".join(entries) + "\n")
    return str(path)


class TestMergeAppliesExclusion:

    def test_listed_accessions_are_dropped_before_the_merge(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCF_2", "GCF_4"])
        rd = RunData()
        selection = _Selection(["GCF_1", "GCF_2", "GCF_3", "GCF_4"])

        merge_wanted_ref_tax(rd, [selection], exclusion_list=excl)

        assert rd.get_wanted_ref_tax_accs() == ["GCF_1", "GCF_3"]

    def test_excluded_count_is_recorded_on_run_data(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCF_2", "GCF_4"])
        rd = RunData()
        merge_wanted_ref_tax(rd, [_Selection(["GCF_1", "GCF_2", "GCF_3", "GCF_4"])],
                             exclusion_list=excl)

        assert rd.wanted_ref_tax_num_excluded == 2

    def test_excluded_count_is_recorded_per_selection(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCF_2"])
        rd = RunData()
        merge_wanted_ref_tax(rd, [_Selection(["GCF_1", "GCF_2", "GCF_3"])],
                             exclusion_list=excl)

        assert rd.wanted_ref_tax_selections[0]["num_excluded"] == 1

    def test_no_exclusion_list_keeps_every_selected_accession(self):
        rd = RunData()
        merge_wanted_ref_tax(rd, [_Selection(["GCF_1", "GCF_2", "GCF_3"])])

        assert rd.get_wanted_ref_tax_accs() == ["GCF_1", "GCF_2", "GCF_3"]
        assert rd.wanted_ref_tax_num_excluded == 0

    def test_accessions_not_in_the_list_are_untouched(self, tmp_path):
        # an exclusion entry that matches nothing removes nothing
        excl = _write(tmp_path / "excl.txt", ["GCF_999"])
        rd = RunData()
        merge_wanted_ref_tax(rd, [_Selection(["GCF_1", "GCF_2"])], exclusion_list=excl)

        assert rd.get_wanted_ref_tax_accs() == ["GCF_1", "GCF_2"]
        assert rd.wanted_ref_tax_num_excluded == 0

    def test_exclusion_spans_multiple_w_selections(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCF_2", "GCF_5"])
        rd = RunData()
        first = _Selection(["GCF_1", "GCF_2"], canonical="Alteromonas")
        second = _Selection(["GCF_4", "GCF_5", "GCF_6"], canonical="Shewanella")

        merge_wanted_ref_tax(rd, [first, second], exclusion_list=excl)

        assert rd.get_wanted_ref_tax_accs() == ["GCF_1", "GCF_4", "GCF_6"]
        assert rd.wanted_ref_tax_num_excluded == 2


class TestMatchingIsOnAccessionCore:
    """
    Exclusion matches on the accession core: GCA_/GCF_ and the version suffix are
    normalised away, so either member of a pair (at any version) excludes the
    assembly however `-w` grabbed it. Real-format accessions are used here so
    accession_core resolves them to a genuine numeric core.
    """

    def test_a_listed_gca_excludes_the_gcf_w_pulled(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCA_000005845.2"])
        rd = RunData()
        # `-w` returned the RefSeq (GCF) member, different version
        merge_wanted_ref_tax(rd, [_Selection(["GCF_000005845.1", "GCF_111111111.1"])],
                             exclusion_list=excl)

        assert rd.get_wanted_ref_tax_accs() == ["GCF_111111111.1"]
        assert rd.wanted_ref_tax_num_excluded == 1

    def test_a_listed_gcf_excludes_the_gca_w_pulled(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCF_000005845.1"])
        rd = RunData()
        merge_wanted_ref_tax(rd, [_Selection(["GCA_000005845.2", "GCA_222222222.1"])],
                             exclusion_list=excl)

        assert rd.get_wanted_ref_tax_accs() == ["GCA_222222222.1"]
        assert rd.wanted_ref_tax_num_excluded == 1

    def test_a_different_version_of_the_same_accession_is_excluded(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCF_000005845.1"])
        rd = RunData()
        merge_wanted_ref_tax(rd, [_Selection(["GCF_000005845.3"])], exclusion_list=excl)

        assert rd.get_wanted_ref_tax_accs() == []
        assert rd.wanted_ref_tax_num_excluded == 1

    def test_a_gtdb_prefixed_entry_normalises_too(self, tmp_path):
        # GTDB's RS_/GB_ prefixes are stripped by accession_core as well
        excl = _write(tmp_path / "excl.txt", ["RS_GCF_000005845.2"])
        rd = RunData()
        merge_wanted_ref_tax(rd, [_Selection(["GCA_000005845.1"])], exclusion_list=excl)

        assert rd.get_wanted_ref_tax_accs() == []
        assert rd.wanted_ref_tax_num_excluded == 1

    def test_a_different_assembly_is_left_alone(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCA_000005845.2"])
        rd = RunData()
        merge_wanted_ref_tax(rd, [_Selection(["GCF_000999999.1"])], exclusion_list=excl)

        assert rd.get_wanted_ref_tax_accs() == ["GCF_000999999.1"]
        assert rd.wanted_ref_tax_num_excluded == 0

    def test_junk_exclusion_entries_never_match(self, tmp_path):
        # entries with no real core ('na', blanks, non-GCA/GCF strings) are dropped
        # from the set, so they can't collide with a `-w` accession that also has an
        # empty core
        excl = _write(tmp_path / "excl.txt", ["na", "not-an-accession"])
        rd = RunData()
        merge_wanted_ref_tax(rd, [_Selection(["some-weird-id", "another-weird-id"])],
                             exclusion_list=excl)

        assert rd.get_wanted_ref_tax_accs() == ["some-weird-id", "another-weird-id"]
        assert rd.wanted_ref_tax_num_excluded == 0


class TestExclusionOnlyTouchesWantedRefTax:

    def test_user_provided_accessions_are_never_excluded(self, tmp_path):
        # a `-a` accession that also appears in the exclusion list stays in the run
        excl = _write(tmp_path / "excl.txt", ["GCF_1"])
        rd = RunData()
        rd.ncbi_accs.append(GenomeData.from_acc("GCF_1"))  # user-provided via -a
        rd.update_all_input_genomes()

        merge_wanted_ref_tax(rd, [_Selection(["GCF_2", "GCF_3"])], exclusion_list=excl)

        assert rd.get_user_provided_ncbi_accs() == ["GCF_1"]
        assert rd.get_wanted_ref_tax_accs() == ["GCF_2", "GCF_3"]


class TestResumeLeavesExclusionCountAlone:

    def test_a_resume_with_no_selections_preserves_the_recorded_count(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCF_2"])
        rd = RunData()
        rd.wanted_ref_tax_num_excluded = 7  # as if a prior run recorded it

        # a resume passes [] because `-w` was already resolved
        merge_wanted_ref_tax(rd, [], exclusion_list=excl)

        assert rd.wanted_ref_tax_num_excluded == 7


class TestPreflightValidation:

    def test_exclusion_list_requires_w(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCF_1"])
        args = _Args(wanted_ref_tax=None, exclusion_list=excl)

        with pytest.raises(SystemExit):
            check_input_genome_files(args)

    def test_a_missing_exclusion_file_exits(self, tmp_path):
        args = _Args(wanted_ref_tax=["Alteromonas"],
                     exclusion_list=str(tmp_path / "nope.txt"))

        with pytest.raises(SystemExit):
            check_input_genome_files(args)

    def test_whitespace_in_the_exclusion_file_exits(self, tmp_path):
        bad = _write(tmp_path / "excl.txt", ["GCF_1 GCF_2"])
        args = _Args(wanted_ref_tax=["Alteromonas"], exclusion_list=bad)

        with pytest.raises(SystemExit):
            check_input_genome_files(args)

    def test_a_clean_exclusion_file_with_w_passes(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCF_1", "GCF_2"])
        args = _Args(wanted_ref_tax=["Alteromonas"], exclusion_list=excl)

        returned = check_input_genome_files(args)

        assert returned.exclusion_list == excl


class TestFingerprint:

    def test_exclusion_list_contents_are_fingerprinted(self, tmp_path):
        excl_a = _write(tmp_path / "a.txt", ["GCF_1"])
        excl_b = _write(tmp_path / "b.txt", ["GCF_1", "GCF_2"])

        base = dict(wanted_ref_tax=["Alteromonas"], ncbi_accessions=None,
                    genbank_files=None, fasta_files=None, amino_acid_files=None,
                    mapping_file=None, target_pfams_file=None, target_kos_file=None,
                    hmm="Universal", target_rank=None, derep_rank="auto",
                    source="gtdb", add_gtdb_tax=False, add_ncbi_tax=False,
                    lineage="", seq_length_cutoff=0.2, gene_representation_cutoff=0.1,
                    genome_hits_cutoff=0.5, best_hit_mode=False, no_super5=False,
                    no_tree=False, tree_program="FastTreeMP", nucleotide_mode=False,
                    keep_gene_alignments=False)

        fp_a = build_fingerprint(_Args(exclusion_list=excl_a, **base))
        fp_b = build_fingerprint(_Args(exclusion_list=excl_b, **base))

        assert fp_a["exclusion_list_sha256"] != fp_b["exclusion_list_sha256"]


class TestReportingMath:

    def test_excluded_line_is_shown_and_kept_out_of_already_counted(self):
        # 10 selected, 2 excluded, 6 added -> 2 were already counted (10 - 6 - 2)
        selection = {
            "taxon": "Alteromonas", "rank": "genus", "derep_rank": "species",
            "num_selected": 10, "num_added": 6, "num_excluded": 2,
        }
        lines = _wanted_ref_tax_detail_lines(selection, indent="")

        assert any("2 selected genomes dropped by the --exclusion-list" in ln
                   for ln in lines)
        assert any("2 more were selected, but were already counted" in ln
                   for ln in lines)

    def test_no_excluded_line_when_none_were_excluded(self):
        selection = {
            "taxon": "Alteromonas", "rank": "genus", "derep_rank": "species",
            "num_selected": 8, "num_added": 8, "num_excluded": 0,
        }
        lines = _wanted_ref_tax_detail_lines(selection, indent="")

        assert not any("--exclusion-list" in ln for ln in lines)
