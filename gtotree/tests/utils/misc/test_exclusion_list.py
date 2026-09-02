"""
`--exclusion-list` keeps named accessions out of what `-w` pulls in.

The exclusion is applied INSIDE the selection core, against the candidate pool, before
dereplication and before any best-per-group pick. That ordering is the whole point:
excluding a species' chosen representative promotes the next-best genome in that
species rather than dropping the species from the run. Post-filtering a dereplicated
set would do the latter, which is what this used to do and deliberately no longer does.

merge_wanted_ref_tax() therefore does NOT filter; it only reads back the count the core
recorded. The file itself is still validated in preflight exactly like an `-a`
accessions file, and is still only meaningful alongside `-w`.
"""

import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore
import pytest  # type: ignore

from gtotree.utils.misc.general import RunData, GenomeData
from gtotree.utils.misc.preflight_checks import (merge_wanted_ref_tax,
                                                 check_input_genome_files,
                                                 build_fingerprint)
from gtotree.utils.misc.messaging import _wanted_ref_tax_detail_lines
from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_derep import select_ref_genomes
from gtotree.utils.taxonomy.exclusion_list import (load_exclusion_cores,
                                                   filter_rows_by_exclusion,
                                                   exclusion_warning)


class _Selection:
    """Stand-in for a RefGenomeSelection as merge/record consume it."""
    def __init__(self, accessions, canonical="Alteromonas",
                 resolved_rank="genus", effective_derep_rank="species",
                 warnings=(), num_excluded=0):
        self.accessions = list(accessions)
        self.canonical = canonical
        self.resolved_rank = resolved_rank
        self.effective_derep_rank = effective_derep_rank
        self.warnings = list(warnings)
        self.num_excluded = num_excluded


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


# ---------------------------------------------------------------------------
# a tiny GTDB-shaped asset: one genus, two species, two genomes each
# ---------------------------------------------------------------------------

_EXTRA_COLS = ["gtdb_representative", "ncbi_refseq_category",
               "checkm2_completeness", "checkm2_contamination",
               "genome_size", "contig_count"]


def _rec(acc, species, comp):
    d = {"ncbi_genbank_assembly_accession": acc,
         "gtdb_representative": "t", "ncbi_refseq_category": "",
         "checkm2_completeness": comp, "checkm2_contamination": "0.5",
         "genome_size": "4000000", "contig_count": "1"}
    lineage = ("Bacteria", "Testophyla", "ClassA", "OrdA", "FamA", "GenA", species)
    for i, r in enumerate(RANKS):
        d[r] = lineage[i]
    return d


# GenA sp1's best genome is ...01 (99.0); its runner-up is ...02 (80.0).
# GenA sp2's best genome is ...03 (99.0); its runner-up is ...04 (70.0).
_RECORDS = [
    _rec("GCA_000000001.1", "GenA sp1", "99.0"),
    _rec("GCA_000000002.1", "GenA sp1", "80.0"),
    _rec("GCA_000000003.1", "GenA sp2", "99.0"),
    _rec("GCA_000000004.1", "GenA sp2", "70.0"),
]


@pytest.fixture
def asset(tmp_path):
    keys = ["ncbi_genbank_assembly_accession"] + list(RANKS) + _EXTRA_COLS
    cols = {k: [rec[k] for rec in _RECORDS] for k in keys}
    path = tmp_path / "gtdb.parquet"
    pq.write_table(
        pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}),
        str(path))
    return str(path)


def _select(asset, exclude=None, derep_rank="species"):
    return select_ref_genomes(asset, "gtdb", "GenA", derep_rank=derep_rank,
                              reps_only=False, exclude_cores=exclude)


class TestExclusionHappensBeforeDereplication:
    """
    The behaviour this feature exists for: an excluded genome leaves the candidate
    pool, so its dereplication group still contributes a genome.
    """

    def test_baseline_picks_the_best_genome_in_each_group(self, asset):
        selection = _select(asset)

        assert selection.accessions == ["GCA_000000001.1", "GCA_000000003.1"]
        assert selection.num_excluded == 0

    def test_excluding_a_groups_best_promotes_its_runner_up(self, asset):
        selection = _select(asset, exclude={"000000001"})

        # sp1 is still represented, by its second-best genome rather than not at all
        assert selection.accessions == ["GCA_000000002.1", "GCA_000000003.1"]

    def test_the_group_is_only_lost_when_every_member_is_excluded(self, asset):
        selection = _select(asset, exclude={"000000001", "000000002"})

        assert selection.accessions == ["GCA_000000003.1"]

    def test_excluded_count_is_candidates_removed_not_genomes_lost(self, asset):
        # two candidates removed, but the selection only shrinks by one group
        selection = _select(asset, exclude={"000000001", "000000002"})

        assert selection.num_excluded == 2

    def test_the_advisory_is_raised_on_the_selection(self, asset):
        selection = _select(asset, exclude={"000000001"})

        assert exclusion_warning(1) in selection.warnings

    def test_derep_off_simply_drops_the_listed_genomes(self, asset):
        selection = _select(asset, exclude={"000000001"}, derep_rank="off")

        assert selection.accessions == ["GCA_000000002.1", "GCA_000000003.1",
                                        "GCA_000000004.1"]
        assert selection.num_excluded == 1

    def test_no_exclusion_keeps_everything(self, asset):
        selection = _select(asset, exclude=None)

        assert len(selection.accessions) == 2
        assert selection.num_excluded == 0
        assert not any("--exclusion-list" in w for w in selection.warnings)


class TestLoadingTheList:

    def test_entries_become_accession_cores(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCF_000005845.2"])

        assert load_exclusion_cores(excl) == {"000005845"}

    def test_no_path_is_an_empty_set(self):
        assert load_exclusion_cores(None) == set()

    def test_blanks_and_comments_are_ignored(self, tmp_path):
        excl = _write(tmp_path / "excl.txt",
                      ["# a comment", "", "GCF_000005845.2", "   "])

        assert load_exclusion_cores(excl) == {"000005845"}

    def test_entries_naming_the_same_assembly_collapse(self, tmp_path):
        excl = _write(tmp_path / "excl.txt",
                      ["GCF_000005845.2", "GCA_000005845.1", "RS_GCF_000005845.3"])

        assert load_exclusion_cores(excl) == {"000005845"}

    def test_junk_entries_are_dropped_rather_than_kept_as_empty_cores(self, tmp_path):
        # otherwise a junk line would match every candidate that also fails to parse
        excl = _write(tmp_path / "excl.txt", ["na", "not-an-accession"])

        assert load_exclusion_cores(excl) == set()


class TestMatchingIsOnAccessionCore:
    """
    GCA_/GCF_ pairing, GTDB's RS_/GB_ prefixes, and version suffixes all normalise
    away, so either member of a pair at any version excludes the assembly however the
    asset happens to spell it.
    """

    def _rows(self):
        return [{"acc": "GCF_000005845.1"}, {"acc": "GCA_111111111.1"}]

    def test_a_listed_gca_excludes_the_gcf_in_the_pool(self, tmp_path):
        cores = load_exclusion_cores(
            _write(tmp_path / "excl.txt", ["GCA_000005845.2"]))

        kept, n = filter_rows_by_exclusion(self._rows(), "acc", cores)

        assert [r["acc"] for r in kept] == ["GCA_111111111.1"]
        assert n == 1

    def test_a_gtdb_prefixed_entry_normalises_too(self, tmp_path):
        cores = load_exclusion_cores(
            _write(tmp_path / "excl.txt", ["RS_GCF_000005845.2"]))

        kept, n = filter_rows_by_exclusion(self._rows(), "acc", cores)

        assert [r["acc"] for r in kept] == ["GCA_111111111.1"]
        assert n == 1

    def test_a_different_assembly_is_left_alone(self, tmp_path):
        cores = load_exclusion_cores(
            _write(tmp_path / "excl.txt", ["GCA_000999999.2"]))

        kept, n = filter_rows_by_exclusion(self._rows(), "acc", cores)

        assert len(kept) == 2
        assert n == 0

    def test_junk_rows_are_not_matched_by_junk_entries(self, tmp_path):
        cores = load_exclusion_cores(_write(tmp_path / "excl.txt", ["na"]))
        rows = [{"acc": "some-weird-id"}, {"acc": "another-weird-id"}]

        kept, n = filter_rows_by_exclusion(rows, "acc", cores)

        assert len(kept) == 2
        assert n == 0


class TestExclusionWarningWording:

    def test_none_excluded_produces_no_line(self):
        assert exclusion_warning(0) is None

    def test_the_line_says_candidate_and_before_selection(self):
        line = exclusion_warning(3)

        assert "3 candidate genomes" in line
        assert "--exclusion-list" in line
        assert "before selection" in line

    def test_one_genome_is_singular(self):
        assert "1 candidate genome " in exclusion_warning(1)


class TestMergeReadsBackTheCount:
    """
    merge_wanted_ref_tax no longer filters; the accessions reaching it are already
    exclusion-free, and it just records what the core reported.
    """

    def test_selection_accessions_are_merged_untouched(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCF_2"])
        rd = RunData()

        merge_wanted_ref_tax(rd, [_Selection(["GCF_1", "GCF_3"], num_excluded=1)],
                             exclusion_list=excl)

        assert rd.get_wanted_ref_tax_accs() == ["GCF_1", "GCF_3"]

    def test_excluded_count_is_recorded_on_run_data(self):
        rd = RunData()
        merge_wanted_ref_tax(rd, [_Selection(["GCF_1"], num_excluded=2)])

        assert rd.wanted_ref_tax_num_excluded == 2

    def test_excluded_count_is_recorded_per_selection(self):
        rd = RunData()
        merge_wanted_ref_tax(rd, [_Selection(["GCF_1"], num_excluded=1)])

        assert rd.wanted_ref_tax_selections[0]["num_excluded"] == 1

    def test_counts_sum_across_multiple_w_selections(self):
        rd = RunData()
        first = _Selection(["GCF_1"], canonical="Alteromonas", num_excluded=1)
        second = _Selection(["GCF_4"], canonical="Shewanella", num_excluded=2)

        merge_wanted_ref_tax(rd, [first, second])

        assert rd.wanted_ref_tax_num_excluded == 3

    def test_no_exclusion_records_zero(self):
        rd = RunData()
        merge_wanted_ref_tax(rd, [_Selection(["GCF_1", "GCF_2", "GCF_3"])])

        assert rd.get_wanted_ref_tax_accs() == ["GCF_1", "GCF_2", "GCF_3"]
        assert rd.wanted_ref_tax_num_excluded == 0


class TestExclusionOnlyTouchesWantedRefTax:

    def test_user_provided_accessions_are_never_excluded(self, tmp_path):
        # a `-a` accession that also appears in the exclusion list stays in the run:
        # the list is only ever handed to the taxonomy selection core
        excl = _write(tmp_path / "excl.txt", ["GCF_1"])
        rd = RunData()
        rd.ncbi_accs.append(GenomeData.from_acc("GCF_1"))  # user-provided via -a
        rd.update_all_input_genomes()

        merge_wanted_ref_tax(rd, [_Selection(["GCF_2", "GCF_3"])],
                             exclusion_list=excl)

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

    def test_excluded_line_is_shown_and_counted_as_candidates(self):
        # the exclusion ran before selection, so num_selected (8) is already net of
        # it: 8 selected, 6 added -> 2 were already counted, independent of the 2
        # candidates the list removed upstream
        selection = {
            "taxon": "Alteromonas", "rank": "genus", "derep_rank": "species",
            "num_selected": 8, "num_added": 6, "num_excluded": 2,
        }
        lines = _wanted_ref_tax_detail_lines(selection, indent="")

        assert any("2 candidate genomes removed by the --exclusion-list "
                   "before selection" in ln for ln in lines)
        assert any("2 more were selected, but were already counted" in ln
                   for ln in lines)

    def test_excluded_genomes_are_not_double_counted_as_overlap(self):
        # everything selected was added, so there is no "already counted" line even
        # though the list removed candidates upstream
        selection = {
            "taxon": "Alteromonas", "rank": "genus", "derep_rank": "species",
            "num_selected": 6, "num_added": 6, "num_excluded": 4,
        }
        lines = _wanted_ref_tax_detail_lines(selection, indent="")

        assert not any("already counted" in ln for ln in lines)

    def test_no_excluded_line_when_none_were_excluded(self):
        selection = {
            "taxon": "Alteromonas", "rank": "genus", "derep_rank": "species",
            "num_selected": 8, "num_added": 8, "num_excluded": 0,
        }
        lines = _wanted_ref_tax_detail_lines(selection, indent="")

        assert not any("--exclusion-list" in ln for ln in lines)
