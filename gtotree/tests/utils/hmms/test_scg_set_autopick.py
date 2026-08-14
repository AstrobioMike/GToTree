"""
Auto-selecting an SCG-set from --wanted-ref-tax (-w) when -H was left off.

The pre-built sets are named for GTDB taxa, so the whole job is: work out where the
request sits in the GTDB taxonomy, then take the finest set covering it.
"""

import pandas as pd  # type: ignore
import pytest  # type: ignore

from gtotree.utils.hmms import scg_hmms_setup as S


# a stand-in for hmm-sources-and-info.tsv holding one of each shape that matters: the
# two domain sets, a phylum, a phylum with GTDB's polyphyly suffix, a class nested
# inside a phylum that also has a set, and the two rankless sets
SUMMARY_ROWS = [
    ("Bacteria.hmm", "bacteria", "domain"),
    ("Archaea.hmm", "archaea", "domain"),
    ("Pseudomonadota.hmm", "pseudomonadota", "phylum"),
    ("Thermoproteota.hmm", "thermoproteota", "phylum"),
    ("Bacillota_I.hmm", "bacillota_i", "phylum"),
    ("Alphaproteobacteria.hmm", "alphaproteobacteria", "class"),
    ("Bacteria-and-Archaea.hmm", "bacteria-and-archaea", "multi-domain"),
    ("Universal-Hug-et-al.hmm", "bacteria,archaea,eukarya", "universal"),
]


@pytest.fixture(autouse=True)
def summary_table(monkeypatch):
    df = pd.DataFrame(SUMMARY_ROWS, columns=["file", "target_taxa", "rank"])
    monkeypatch.setattr(S, "read_in_hmm_summary_table", lambda: df)


class FakeSelection:
    """Enough of a RefGenomeSelection for the auto-pick to work on."""

    def __init__(self, canonical, resolved_rank, rows=None, accessions=None):
        self.canonical = canonical
        self.resolved_rank = resolved_rank
        self.rows = rows or []
        self.accessions = accessions or []


def lineage(**ranks):
    return ranks


################################################################################
# indexing the summary table
################################################################################

def test_only_sets_that_sit_at_a_real_rank_are_indexed():
    index = S.packaged_sets_by_taxon()

    assert index[("domain", "bacteria")] == "Bacteria"
    assert index[("phylum", "pseudomonadota")] == "Pseudomonadota"
    assert index[("class", "alphaproteobacteria")] == "Alphaproteobacteria"

    # 'universal' and 'multi-domain' aren't nodes in a lineage, so a walk must never
    # be able to land on them
    assert all(rank in S.RANKS for rank, _ in index)
    assert "Universal-Hug-et-al" not in index.values()
    assert "Bacteria-and-Archaea" not in index.values()


def test_an_unreadable_summary_table_yields_an_empty_index(monkeypatch):
    monkeypatch.setattr(S, "read_in_hmm_summary_table",
                        lambda: (_ for _ in ()).throw(OSError("nope")))
    assert S.packaged_sets_by_taxon() == {}


################################################################################
# walking a lineage
################################################################################

def test_a_taxon_with_its_own_set_matches_at_its_own_rank():
    full = lineage(domain="Bacteria", phylum="Pseudomonadota",
                   **{"class": "Alphaproteobacteria"})
    assert S.pick_set_for_lineage(full, "phylum") == (
        "Pseudomonadota", "phylum", "Pseudomonadota")


def test_a_taxon_without_a_set_climbs_to_its_nearest_ancestor_that_has_one():
    full = lineage(domain="Bacteria", phylum="Pseudomonadota",
                   **{"class": "Gammaproteobacteria"}, order="Enterobacterales")
    assert S.pick_set_for_lineage(full, "order") == (
        "Pseudomonadota", "phylum", "Pseudomonadota")


def test_the_finest_covering_set_wins_over_a_coarser_one():
    full = lineage(domain="Bacteria", phylum="Pseudomonadota",
                   **{"class": "Alphaproteobacteria"}, order="Rhizobiales")
    assert S.pick_set_for_lineage(full, "order")[0] == "Alphaproteobacteria"


def test_ranks_below_the_resolved_rank_are_ignored():
    """
    A row pulled for a phylum-level request still carries genus and species values --
    whichever genome happened to be read. Honouring them would pick a set for an
    arbitrary sliver of the request.
    """
    full = lineage(domain="Bacteria", phylum="Bacteroidota",
                   **{"class": "Alphaproteobacteria"})
    assert S.pick_set_for_lineage(full, "phylum") == ("Bacteria", "domain", "Bacteria")


def test_ranks_missing_from_a_lineage_are_skipped():
    full = lineage(domain="Bacteria", phylum="NA", **{"class": ""})
    assert S.pick_set_for_lineage(full, "class") == ("Bacteria", "domain", "Bacteria")


def test_nothing_matches_when_no_ancestor_has_a_set():
    full = lineage(domain="Eukaryota", phylum="Ascomycota")
    assert S.pick_set_for_lineage(full, "phylum") == (None, None, None)


def test_an_unknown_resolved_rank_matches_nothing():
    assert S.pick_set_for_lineage(lineage(domain="Bacteria"), "kingdom") == (None, None, None)


################################################################################
# collapsing many lineages into one
################################################################################

def test_unanimous_lineages_collapse_to_themselves():
    one = lineage(domain="Bacteria", phylum="Pseudomonadota",
                  **{"class": "Alphaproteobacteria"})
    agreed, deepest = S.consensus_lineage([one, dict(one), dict(one)])

    assert deepest == "class"
    assert agreed["phylum"] == "Pseudomonadota"


def test_agreement_ends_at_the_first_rank_the_genomes_disagree_on():
    lineages = [lineage(domain="Bacteria", phylum="Pseudomonadota",
                        **{"class": "Alphaproteobacteria"}),
                lineage(domain="Bacteria", phylum="Pseudomonadota",
                        **{"class": "Gammaproteobacteria"})]
    agreed, deepest = S.consensus_lineage(lineages)

    assert deepest == "phylum"
    assert "class" not in agreed


def test_deep_agreement_below_a_broken_rank_is_not_used():
    """
    The reason consensus walks coarse-to-fine with an early stop. These genomes split
    at phylum, but the two that are classified to genus happen to share one -- taking
    the deepest unanimous rank would pick a set covering a sliver of the request.
    """
    lineages = [lineage(domain="Bacteria", phylum="Pseudomonadota", genus="Escherichia"),
                lineage(domain="Bacteria", phylum="Bacteroidota", genus="NA"),
                lineage(domain="Bacteria", phylum="Bacillota", genus="Escherichia")]
    agreed, deepest = S.consensus_lineage(lineages)

    assert deepest == "domain"
    assert "genus" not in agreed


def test_a_small_minority_does_not_break_agreement():
    lineages = ([lineage(domain="Bacteria", phylum="Bacillota")] * 19
                + [lineage(domain="Bacteria", phylum="Bacillota_I")])
    agreed, deepest = S.consensus_lineage(lineages)

    assert deepest == "phylum"
    assert agreed["phylum"] == "Bacillota"


def test_genomes_from_different_domains_agree_on_nothing():
    lineages = [lineage(domain="Bacteria", phylum="Pseudomonadota"),
                lineage(domain="Archaea", phylum="Thermoproteota")]
    assert S.consensus_lineage(lineages) == ({}, None)


def test_no_lineages_at_all_agree_on_nothing():
    assert S.consensus_lineage([]) == ({}, None)


################################################################################
# the GTDB path
################################################################################

def test_a_gtdb_taxon_with_its_own_set_gets_it():
    selection = FakeSelection(
        "Pseudomonadota", "phylum",
        rows=[lineage(domain="Bacteria", phylum="Pseudomonadota", genus="Escherichia")])
    picked = S.autopick_scg_set("gtdb", selection)

    assert picked.name == "Pseudomonadota"


def test_a_gtdb_taxon_below_any_set_climbs_to_its_ancestor():
    selection = FakeSelection(
        "Rhizobiales", "order",
        rows=[lineage(domain="Bacteria", phylum="Pseudomonadota",
                      **{"class": "Alphaproteobacteria"}, order="Rhizobiales")])
    picked = S.autopick_scg_set("gtdb", selection)

    assert picked.name == "Alphaproteobacteria"
    assert picked.matched_rank == "class"
    assert "sits within class Alphaproteobacteria" in picked.reason


def test_the_gtdb_path_needs_no_lookup_beyond_the_selection(monkeypatch):
    """The selection's rows already carry the GTDB lineage."""
    import gtotree.utils.gtdb.handle_gtdb_tax_info as G
    monkeypatch.setattr(G, "gtdb_lineages_for_accessions",
                        lambda *a, **k: pytest.fail("should not link on the GTDB path"))

    selection = FakeSelection("Thermoproteota", "phylum",
                              rows=[lineage(domain="Archaea", phylum="Thermoproteota")],
                              accessions=["GCA_000001.1"])
    assert S.autopick_scg_set("gtdb", selection).name == "Thermoproteota"


################################################################################
# the NCBI path
################################################################################

@pytest.fixture
def linked(monkeypatch):
    """Stand in for the accession -> GTDB lineage join."""
    import gtotree.utils.gtdb.handle_gtdb_tax_info as G

    def set_result(lineages):
        monkeypatch.setattr(
            G, "gtdb_lineages_for_accessions",
            lambda accessions, gtdb_path=None: {
                f"GCA_{i}": lin for i, lin in enumerate(lineages)})

    return set_result


def test_an_ncbi_taxon_is_matched_through_its_genomes_gtdb_placement(linked):
    # NCBI's Nitrososphaerota is GTDB's Thermoproteota -- no name in common
    linked([lineage(domain="Archaea", phylum="Thermoproteota")] * 8)
    picked = S.autopick_scg_set("ncbi", FakeSelection("Nitrososphaerota", "phylum",
                                                      accessions=["GCA_x"] * 8))

    assert picked.name == "Thermoproteota"
    assert "fall in phylum Thermoproteota in GTDB" in picked.reason


def test_a_target_gtdb_has_never_heard_of_falls_back_to_universal(linked):
    linked([])
    picked = S.autopick_scg_set("ncbi", FakeSelection("Ascomycota", "phylum",
                                                      accessions=["GCA_x"]))

    assert picked.name == S.UNIVERSAL_SCG_SET
    assert "counterpart in GTDB" in picked.reason


def test_a_target_straddling_the_domains_takes_the_cross_domain_set(linked):
    linked([lineage(domain="Bacteria", phylum="Pseudomonadota"),
            lineage(domain="Archaea", phylum="Thermoproteota")])
    picked = S.autopick_scg_set("ncbi", FakeSelection("Whatever", "phylum",
                                                      accessions=["a", "b"]))

    assert picked.name == S.CROSS_DOMAIN_SCG_SET
    assert "both domains" in picked.reason


def test_genomes_that_agree_only_on_a_domain_get_that_domains_set(linked):
    linked([lineage(domain="Bacteria", phylum="Pseudomonadota"),
            lineage(domain="Bacteria", phylum="Bacteroidota")])
    assert S.autopick_scg_set("ncbi", FakeSelection("Whatever", "phylum",
                                                    accessions=["a", "b"])).name == "Bacteria"


################################################################################
# no answer is a fallback, never an error
################################################################################

def test_no_selection_at_all_falls_back_to_universal():
    assert S.autopick_scg_set("gtdb", None).name == S.UNIVERSAL_SCG_SET


def test_an_empty_summary_table_falls_back_to_universal(monkeypatch):
    monkeypatch.setattr(S, "packaged_sets_by_taxon", lambda: {})
    selection = FakeSelection("Pseudomonadota", "phylum",
                              rows=[lineage(domain="Bacteria", phylum="Pseudomonadota")])

    assert S.autopick_scg_set("gtdb", selection).name == S.UNIVERSAL_SCG_SET


################################################################################
# several -w taxa at once
################################################################################

def test_two_taxa_under_one_phylum_get_that_phylums_set():
    picked = S.autopick_scg_set("gtdb", [
        FakeSelection("Alphaproteobacteria", "class",
                      rows=[lineage(domain="Bacteria", phylum="Pseudomonadota",
                                    **{"class": "Alphaproteobacteria"})]),
        FakeSelection("Gammaproteobacteria", "class",
                      rows=[lineage(domain="Bacteria", phylum="Pseudomonadota",
                                    **{"class": "Gammaproteobacteria"})]),
    ])

    assert picked.name == "Pseudomonadota"
    assert "sit within phylum Pseudomonadota" in picked.reason


def test_two_phyla_in_one_domain_climb_to_that_domains_set():
    picked = S.autopick_scg_set("gtdb", [
        FakeSelection("Pseudomonadota", "phylum",
                      rows=[lineage(domain="Bacteria", phylum="Pseudomonadota")]),
        FakeSelection("Bacillota_I", "phylum",
                      rows=[lineage(domain="Bacteria", phylum="Bacillota_I")]),
    ])

    assert picked.name == "Bacteria"


def test_bacteria_plus_archaea_takes_the_cross_domain_set():
    picked = S.autopick_scg_set("gtdb", [
        FakeSelection("Bacteria", "domain", rows=[lineage(domain="Bacteria")]),
        FakeSelection("Archaea", "domain", rows=[lineage(domain="Archaea")]),
    ])

    assert picked.name == S.CROSS_DOMAIN_SCG_SET
    assert "both domains" in picked.reason


def test_a_broad_taxon_is_not_dragged_below_its_own_rank_by_a_narrow_one():
    """
    `-w Bacteria -w Pseudomonadota` has to cover Bacteria, so the answer is the domain
    set -- a taxon can't constrain ranks below the one it resolved to.
    """
    picked = S.autopick_scg_set("gtdb", [
        FakeSelection("Bacteria", "domain", rows=[lineage(domain="Bacteria")]),
        FakeSelection("Pseudomonadota", "phylum",
                      rows=[lineage(domain="Bacteria", phylum="Pseudomonadota")]),
    ])

    assert picked.name == "Bacteria"


def test_a_single_selection_still_works_unwrapped():
    """One selection, passed bare rather than in a list, as the driver may still do."""
    selection = FakeSelection(
        "Pseudomonadota", "phylum",
        rows=[lineage(domain="Bacteria", phylum="Pseudomonadota")])

    assert S.autopick_scg_set("gtdb", selection).name == "Pseudomonadota"
    assert S.autopick_scg_set("gtdb", [selection]).name == "Pseudomonadota"


def test_an_empty_list_of_selections_falls_back_to_universal():
    assert S.autopick_scg_set("gtdb", []).name == S.UNIVERSAL_SCG_SET


################################################################################
# eukaryotes are the one thing that still takes the universal set
################################################################################

def test_a_eukaryotic_selection_falls_back_to_universal():
    picked = S.autopick_scg_set("ncbi", FakeSelection(
        "Ascomycota", "phylum",
        rows=[lineage(domain="Eukaryota", phylum="Ascomycota")],
        accessions=["GCA_x"]))

    assert picked.name == S.UNIVERSAL_SCG_SET
    assert "outside Bacteria and Archaea" in picked.reason


def test_eukaryotes_alongside_bacteria_still_take_the_universal_set(linked):
    """
    The case the GTDB join would hide: linking drops the eukaryotes, and what's left
    would otherwise look like a clean bacterial pick.
    """
    linked([lineage(domain="Bacteria", phylum="Pseudomonadota")] * 20)
    picked = S.autopick_scg_set("ncbi", [
        FakeSelection("Pseudomonadota", "phylum",
                      rows=[lineage(domain="Bacteria", phylum="Pseudomonadota")],
                      accessions=["GCA_x"] * 20),
        FakeSelection("Ascomycota", "phylum",
                      rows=[lineage(domain="Eukaryota", phylum="Ascomycota")],
                      accessions=["GCA_y"]),
    ])

    assert picked.name == S.UNIVERSAL_SCG_SET
    assert "outside Bacteria and Archaea" in picked.reason


def test_a_gtdb_selection_never_looks_eukaryotic():
    """GTDB classifies only Bacteria and Archaea, so the check is a no-op there."""
    assert S.domains_outside_prebuilt_scope(
        [lineage(domain="Bacteria"), lineage(domain="Archaea")]) == set()
    assert S.domains_outside_prebuilt_scope([lineage(domain="NA"), lineage()]) == set()


################################################################################
# intersecting independently resolved taxa
################################################################################

def test_intersecting_taxa_is_unanimous_not_a_vote():
    """
    A 90% vote would let the big taxon swallow the small one. Nineteen bacterial taxa
    and one archaeal one still have to come out cross-domain.
    """
    resolved = ([(lineage(domain="Bacteria", phylum="Bacillota"), "phylum")] * 19
                + [(lineage(domain="Archaea", phylum="Thermoproteota"), "phylum")])

    assert S.intersect_lineages(resolved) == ({}, None)


def test_intersecting_stops_at_the_coarsest_taxons_rank():
    resolved = [(lineage(domain="Bacteria"), "domain"),
                (lineage(domain="Bacteria", phylum="Pseudomonadota"), "phylum")]
    agreed, deepest = S.intersect_lineages(resolved)

    assert deepest == "domain"
    assert "phylum" not in agreed


def test_intersecting_nothing_agrees_on_nothing():
    assert S.intersect_lineages([]) == ({}, None)


def test_intersecting_an_unknown_rank_agrees_on_nothing():
    assert S.intersect_lineages([(lineage(domain="Bacteria"), "kingdom")]) == ({}, None)


################################################################################
# viral targets bring their own HMMs
################################################################################

def viral_args(hmm=None):
    import argparse
    return argparse.Namespace(hmm=hmm)


def test_a_viral_selection_is_spotted_from_its_own_rows():
    selection = FakeSelection("Caudoviricetes", "class",
                              rows=[lineage(domain="Viruses", phylum="Uroviricota")])

    assert S.viral_selections(selection) == [selection]


@pytest.mark.parametrize("domain", ["Viruses", "viruses", "Viroids", "viral metagenome"])
def test_the_viral_domain_names_that_count(domain):
    selection = FakeSelection("Something", "phylum", rows=[lineage(domain=domain)])

    assert S.viral_selections(selection)


@pytest.mark.parametrize("domain", ["Bacteria", "Archaea", "Eukaryota", "NA", ""])
def test_nothing_cellular_is_mistaken_for_viral(domain):
    selection = FakeSelection("Something", "phylum", rows=[lineage(domain=domain)])

    assert S.viral_selections(selection) == []


def test_a_viral_target_with_no_hmm_given_is_refused():
    selection = FakeSelection("Caudoviricetes", "class",
                              rows=[lineage(domain="Viruses")])

    with pytest.raises(S.ViralTaxonNeedsOwnHMMs) as excinfo:
        S.check_viruses_have_their_own_hmms(viral_args(), selection)

    assert "Caudoviricetes" in str(excinfo.value)
    assert "-H" in str(excinfo.value)


def test_a_viral_target_cannot_take_a_pre_built_set_either():
    """
    Naming a pre-built set is refused just as flatly as naming none: the point is that
    no set built from bacterial and archaeal genomes means anything for a viral tree.
    """
    selection = FakeSelection("Caudoviricetes", "class",
                              rows=[lineage(domain="Viruses")])

    with pytest.raises(S.ViralTaxonNeedsOwnHMMs) as excinfo:
        S.check_viruses_have_their_own_hmms(viral_args(hmm="Universal-Hug-et-al"),
                                            selection)

    assert "Universal-Hug-et-al" in str(excinfo.value)


def test_a_viral_target_with_the_users_own_hmm_file_is_allowed(tmp_path):
    hmm_file = tmp_path / "my-viral-SCGs.hmm"
    hmm_file.write_text("HMMER3/f\n")

    selection = FakeSelection("Caudoviricetes", "class",
                              rows=[lineage(domain="Viruses")])

    S.check_viruses_have_their_own_hmms(viral_args(hmm=str(hmm_file)), selection)


def test_one_viral_taxon_among_several_is_enough_to_refuse():
    bacterial = FakeSelection("Pseudomonadota", "phylum",
                              rows=[lineage(domain="Bacteria", phylum="Pseudomonadota")])
    viral = FakeSelection("Caudoviricetes", "class", rows=[lineage(domain="Viruses")])

    with pytest.raises(S.ViralTaxonNeedsOwnHMMs) as excinfo:
        S.check_viruses_have_their_own_hmms(viral_args(), [bacterial, viral])

    # only the offending taxon is named, not the innocent one alongside it
    assert "Caudoviricetes" in str(excinfo.value)
    assert "Pseudomonadota" not in str(excinfo.value)


def test_nothing_viral_means_nothing_to_check():
    selection = FakeSelection("Pseudomonadota", "phylum",
                              rows=[lineage(domain="Bacteria", phylum="Pseudomonadota")])

    S.check_viruses_have_their_own_hmms(viral_args(), selection)
    S.check_viruses_have_their_own_hmms(viral_args(), None)


def test_a_set_named_for_the_taxon_itself_gives_no_reason():
    """
    The banner prints the taxon and the set it picked on the line above, so a reason
    saying 'Bacteria' got the 'Bacteria' set would only repeat it. Empty means the
    banner falls back to a bare '(auto-selected)'.
    """
    selection = FakeSelection(
        "Pseudomonadota", "phylum",
        rows=[lineage(domain="Bacteria", phylum="Pseudomonadota")])
    picked = S.autopick_scg_set("gtdb", selection)

    assert picked.name == "Pseudomonadota"
    assert picked.reason == ""
