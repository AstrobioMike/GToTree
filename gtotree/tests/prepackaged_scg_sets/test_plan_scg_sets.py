import importlib.util
from pathlib import Path
import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore
import pytest  # type: ignore
from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_derep import select_ref_genomes


# anchored to this file rather than counted up with parents[N], which breaks whenever
# a test file moves to a different depth
REPO_ROOT = Path(__file__).resolve().parents[3]
PLANNER_PATH = REPO_ROOT / "hmm_sets" / "creation" / "plan-scg-sets.py"


@pytest.fixture(scope="module")
def planner():
    """Load the dev script by path -- it isn't importable as a package module."""
    if not PLANNER_PATH.is_file():
        pytest.skip(f"planner script not found at {PLANNER_PATH}")
    spec = importlib.util.spec_from_file_location("plan_scg_sets", PLANNER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GTDB_COLS = ["ncbi_genbank_assembly_accession", "gtdb_representative",
             "ncbi_refseq_category", "checkm2_completeness", "checkm2_contamination",
             "genome_size", "contig_count"]


@pytest.fixture(scope="module")
def mock_gtdb(tmp_path_factory):
    """
    A small GTDB-shaped parquet covering the cases the planner has to get right:

      Bigphylum   -- several genera, a mix of quality, one genus (GenusD) that exists
                     only below the floor and so must not be counted
      Smallphylum -- few genera but multiple species per genus, so genus-level derep
                     under-delivers and species-level derep is the right answer
      a non-representative row, which the reps-only default must exclude
      a row with no genus assigned, which grouping must skip
    """
    path = tmp_path_factory.mktemp("planner") / "gtdb-data.parquet"
    records = []

    def add(acc, lineage, comp, cont, rep="t"):
        records.append((acc, lineage, comp, cont, rep))

    def lin(phylum, cls, genus, species):
        return ("Bacteria", phylum, cls, f"{cls}ales", f"{cls}aceae", genus, species)

    # Bigphylum: GenusA/B/C clear the floor, GenusD does not
    add("GCA_001.1", lin("Bigphylum", "Bigia", "GenusA", "GenusA sp1"), "99.0", "0.5")
    add("GCA_002.1", lin("Bigphylum", "Bigia", "GenusA", "GenusA sp2"), "97.0", "1.0")
    add("GCA_003.1", lin("Bigphylum", "Bigia", "GenusB", "GenusB sp1"), "95.0", "2.0")
    add("GCA_004.1", lin("Bigphylum", "Bigib", "GenusC", "GenusC sp1"), "91.0", "4.0")
    add("GCA_005.1", lin("Bigphylum", "Bigib", "GenusD", "GenusD sp1"), "60.0", "1.0")
    add("GCA_006.1", lin("Bigphylum", "Bigib", "GenusD", "GenusD sp2"), "70.0", "9.0")
    # not a representative: excluded by the reps-only default even though it's pristine
    add("GCA_007.1", lin("Bigphylum", "Bigia", "GenusE", "GenusE sp1"), "99.9", "0.0",
        rep="f")
    # no genus assigned: must be skipped by genus-level grouping
    add("GCA_008.1", lin("Bigphylum", "Bigia", "", "unclassified sp1"), "99.0", "0.5")

    # Smallphylum: 2 genera but 5 species, the shape that forces a species-level derep
    add("GCA_101.1", lin("Smallphylum", "Smallia", "GenusS", "GenusS sp1"), "99.0", "0.5")
    add("GCA_102.1", lin("Smallphylum", "Smallia", "GenusS", "GenusS sp2"), "98.0", "0.5")
    add("GCA_103.1", lin("Smallphylum", "Smallia", "GenusS", "GenusS sp3"), "97.0", "0.5")
    add("GCA_104.1", lin("Smallphylum", "Smallia", "GenusT", "GenusT sp1"), "96.0", "0.5")
    add("GCA_105.1", lin("Smallphylum", "Smallia", "GenusT", "GenusT sp2"), "95.0", "0.5")

    cols = {c: [] for c in GTDB_COLS}
    for r in RANKS:
        cols[r] = []
    for acc, lineage, comp, cont, rep in records:
        cols["ncbi_genbank_assembly_accession"].append(acc)
        cols["gtdb_representative"].append(rep)
        cols["ncbi_refseq_category"].append("na")
        cols["checkm2_completeness"].append(comp)
        cols["checkm2_contamination"].append(cont)
        cols["genome_size"].append("4000000")
        cols["contig_count"].append("50")
        for i, rank in enumerate(RANKS):
            cols[rank].append(lineage[i])

    pq.write_table(
        pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}),
        str(path))
    return str(path)


@pytest.fixture(scope="module")
def stats(planner, mock_gtdb):
    return planner.gather_stats(mock_gtdb, 90.0, 5.0, ranks=("phylum", "class"))


################################################################################
# counting
################################################################################

def test_quality_floor_removes_a_genus_from_the_counts(stats):
    """GenusD has only sub-floor genomes, so it must not be counted."""
    assert stats[("phylum", "Bigphylum")].n("genus") == 3


def test_non_representatives_are_not_counted(stats):
    """GenusE's only genome is pristine but isn't a GTDB representative, and the
    selector's reps-only default won't see it -- so neither should the planner."""
    genera = stats[("phylum", "Bigphylum")].groups["genus"]
    assert "GenusE" not in genera


def test_unassigned_rank_values_are_skipped(stats):
    """derep() skips rows with no value at the grouping rank, so counting them would
    overpromise by one."""
    assert "" not in stats[("phylum", "Bigphylum")].groups["genus"]


def test_total_genome_count_ignores_the_floor(stats):
    """The genome count is the sequencing-effort signal, so it's deliberately taken
    before any filtering -- unlike the genus count."""
    assert stats[("phylum", "Bigphylum")].genomes == 8


def test_class_level_stats_are_gathered_independently(stats):
    # Bigia holds GenusA and GenusB (GenusE isn't a rep, and one row has no genus)
    assert stats[("class", "Bigia")].n("genus") == 2
    # Bigib holds GenusC and GenusD, but GenusD is entirely below the floor
    assert stats[("class", "Bigib")].n("genus") == 1


################################################################################
# the drift guard: predictions vs. the real selector
################################################################################

@pytest.mark.parametrize("rank,taxon,derep", [
    ("phylum", "Bigphylum", "genus"),
    ("phylum", "Bigphylum", "species"),
    ("phylum", "Smallphylum", "genus"),
    ("phylum", "Smallphylum", "species"),
    ("class", "Bigia", "genus"),
])
def test_predicted_counts_match_what_the_selector_returns(stats, mock_gtdb,
                                                          rank, taxon, derep):
    predicted = stats[(rank, taxon)].n(derep)
    selection = select_ref_genomes(mock_gtdb, "gtdb", taxon, target_rank=rank,
                                   derep_rank=derep, min_completeness=90.0,
                                   max_contamination=5.0)
    assert predicted == len(selection.accessions)


################################################################################
# derep rank choice
################################################################################

def test_derep_steps_finer_when_the_preferred_rank_is_too_small(planner, stats):
    """Smallphylum has 2 genera but 5 species; a min of 4 should force species."""
    rank, count, flags = planner.choose_derep_rank(
        stats[("phylum", "Smallphylum")], "genus", 4, None)
    assert rank == "species" and count == 5
    assert any("stepped" in f for f in flags)


def test_derep_stays_at_the_preferred_rank_when_it_suffices(planner, stats):
    rank, count, flags = planner.choose_derep_rank(
        stats[("phylum", "Smallphylum")], "genus", 2, None)
    assert rank == "genus" and count == 2 and not flags


def test_derep_never_steps_coarser_on_its_own(planner, stats):
    """A large input is flagged for a human rather than silently coarsened, because
    coarsening past genus throws away real breadth."""
    rank, count, flags = planner.choose_derep_rank(
        stats[("phylum", "Bigphylum")], "genus", 1, wants_review_above=1)
    assert rank == "genus"
    assert any("coarser" in f for f in flags)


def test_unreachable_minimum_is_flagged_not_hidden(planner, stats):
    """If even species-level derep can't reach the minimum, the plan says so rather
    than quietly promising a weak set."""
    rank, count, flags = planner.choose_derep_rank(
        stats[("phylum", "Smallphylum")], "genus", 500, None)
    assert rank == "species"
    assert any("weak" in f for f in flags)


################################################################################
# plan assembly
################################################################################

def _config(**overrides):
    config = {
        "quality": {"min_completeness": 90.0, "max_contamination": 5.0},
        "criteria": {"rank": "phylum", "min_genera": 3, "min_genomes": 1000},
        "derep": {"preferred_rank": "genus", "min_input_genomes": 2},
        "include": [],
        "exclude": [],
    }
    config.update(overrides)
    return config


def test_either_criterion_selects_a_phylum(planner, stats):
    """The genera and genome rules are an OR: each catches clades the other misses."""
    plan, _ = planner.build_plan(_config(), stats)
    names = {p.name for p in plan}
    assert "Bigphylum" in names        # 3 genera, clears min_genera
    assert "Smallphylum" not in names  # 2 genera, 5 genomes: clears neither


def test_genome_count_alone_can_select(planner, stats):
    plan, _ = planner.build_plan(
        _config(criteria={"rank": "phylum", "min_genera": 999, "min_genomes": 5}), stats)
    assert {p.name for p in plan} == {"Bigphylum", "Smallphylum"}


def test_include_adds_a_taxon_below_the_bars(planner, stats):
    plan, _ = planner.build_plan(
        _config(include=[{"rank": "phylum", "name": "Smallphylum"}]), stats)
    entry = next(p for p in plan if p.name == "Smallphylum")
    assert entry.reason == "include"


def test_include_reaches_ranks_automatic_selection_never_does(planner, stats):
    """Automatic selection runs at phylum only, so a class can only arrive via include."""
    plan, _ = planner.build_plan(
        _config(include=[{"rank": "class", "name": "Bigia"}]), stats)
    assert any(p.rank == "class" and p.name == "Bigia" for p in plan)


def test_include_on_an_auto_selected_taxon_annotates_rather_than_duplicates(planner, stats):
    plan, _ = planner.build_plan(
        _config(include=[{"rank": "phylum", "name": "Bigphylum"}]),
        stats)
    matches = [p for p in plan if p.name == "Bigphylum"]
    assert len(matches) == 1
    assert "include" in matches[0].reason


def test_exclude_beats_both_criteria_and_include(planner, stats):
    plan, _ = planner.build_plan(
        _config(include=[{"rank": "phylum", "name": "Bigphylum"}],
                exclude=[{"rank": "phylum", "name": "Bigphylum"}]), stats)
    assert not any(p.name == "Bigphylum" for p in plan)


def test_configured_taxon_absent_from_gtdb_is_reported(planner, stats):
    """A renamed or reclassified taxon should surface loudly, not vanish from the plan."""
    plan, missing = planner.build_plan(
        _config(include=[{"rank": "phylum", "name": "Gonephylum"}]), stats)
    assert missing and "Gonephylum" in missing[0]
    assert not any(p.name == "Gonephylum" for p in plan)


################################################################################
# domain-spanning sets
################################################################################

def _domain_config(**overrides):
    config = {
        "quality": {"min_completeness": 90.0, "max_contamination": 5.0},
        "domain_set": [
            {"name": "Bacteria", "taxa": ["Bacteria"], "derep_rank": "family"},
            {"name": "Archaea", "taxa": ["Archaea"], "derep_rank": "family"},
            {"name": "Bacteria-and-Archaea", "taxa": ["Bacteria", "Archaea"],
             "derep_rank": "family"},
        ],
    }
    config.update(overrides)
    return config


def test_domain_sets_are_planned_from_their_block(planner, mock_gtdb):
    plan = planner.plan_domain_sets(_domain_config(), mock_gtdb,
                                    {"min_completeness": 90.0, "max_contamination": 5.0},
                                    None)
    names = {p.name for p in plan}
    assert names == {"Bacteria", "Archaea", "Bacteria-and-Archaea"}


def test_single_domain_set_counts_families_in_that_domain(planner, mock_gtdb):
    """The mock has only Bacteria, whose two phyla contribute distinct families."""
    plan = planner.plan_domain_sets(_domain_config(), mock_gtdb,
                                    {"min_completeness": 90.0, "max_contamination": 5.0},
                                    None)
    bacteria = next(p for p in plan if p.name == "Bacteria")
    # Bigia/Bigib/Smallia families, minus GenusD's (floored out) -- Bigiaaceae,
    # Bigibaceae, Smalliaaceae
    assert bacteria.n_input == 3
    assert bacteria.derep_rank == "family"


def test_combined_set_sums_across_its_domains(planner, mock_gtdb):
    """The mock has no Archaea, so Bacteria-and-Archaea equals Bacteria alone here --
    the point is that the count is a sum over the named domains, not a separate query."""
    plan = planner.plan_domain_sets(_domain_config(), mock_gtdb,
                                    {"min_completeness": 90.0, "max_contamination": 5.0},
                                    None)
    bacteria = next(p for p in plan if p.name == "Bacteria")
    combined = next(p for p in plan if p.name == "Bacteria-and-Archaea")
    assert combined.n_input == bacteria.n_input   # Archaea contributes 0 in the mock
    assert combined.taxa == ["Bacteria", "Archaea"]
    assert combined.rank == "multi-domain"


def test_domain_set_predicted_count_matches_selector_union(planner, mock_gtdb):
    """The drift guard, for the combined set: the planner's sum-over-domains must equal
    the union the build produces from repeated -w."""
    from gtotree.utils.taxonomy.tax_derep import select_ref_genomes
    plan = planner.plan_domain_sets(_domain_config(), mock_gtdb,
                                    {"min_completeness": 90.0, "max_contamination": 5.0},
                                    None)
    bacteria = next(p for p in plan if p.name == "Bacteria")
    selection = select_ref_genomes(mock_gtdb, "gtdb", "Bacteria", target_rank="domain",
                                   derep_rank="family", min_completeness=90.0,
                                   max_contamination=5.0)
    assert bacteria.n_input == len(selection.accessions)


def test_percent_single_copy_override_is_carried(planner, mock_gtdb):
    plan = planner.plan_domain_sets(
        _domain_config(**{"domain_set": [
            {"name": "Bacteria", "taxa": ["Bacteria"], "derep_rank": "family",
             "percent_single_copy": 75}]}),
        mock_gtdb, {"min_completeness": 90.0, "max_contamination": 5.0}, None)
    assert plan[0].percent_single_copy == 75


def test_no_domain_block_plans_nothing(planner, mock_gtdb):
    plan = planner.plan_domain_sets({}, mock_gtdb,
                                    {"min_completeness": 90.0, "max_contamination": 5.0},
                                    None)
    assert plan == []
