"""
The `-w` drivers must hand `--gtdb-section` to the selection core, and must hand the
empty-result messaging THEIR OWN flag names.

`resolve_input_genomes` is shared by the main GToTree driver, `gtt gen-scg-hmms` and
`gtt search-annotations`, so one set of assertions covers all three surfaces.

The second half of this module is the trap the flag walks into. Before it existed the
drivers passed `reps_only=None`, so `empty_selection_message` received
`reps_only_requested=False` and took its "this is the source's default" branch. Start
passing `reps_only=True` and the naive thing -- letting `reps_only_requested` follow
the effective value -- flips it to the "the user asked for this" branch, which names
`--representatives-only`: a flag that lives on `gtt dl-ncbi-assemblies` and on no `-w`
driver at all. An empty pull would tell people to loosen something they cannot spell.
"""

import pytest  # type: ignore

from gtotree.utils.misc.general import RunData, resolve_input_genomes
from gtotree.utils.taxonomy.tax_derep import PoolAttrition, RefGenomeSelection


class _Args:
    def __init__(self, **kw):
        self.wanted_ref_tax = ["Testophyla"]
        self.source = "gtdb"
        self.gtdb_section = None
        self.ncbi_section = "both"
        self.target_rank = None
        self.target_domain = None
        self.derep_rank = "auto"
        self.exclusion_list = None
        self.ncbi_accessions = None
        self.genbank_files = None
        self.fasta_files = None
        self.amino_acid_files = None
        self.__dict__.update(kw)


class _Error(Exception):
    pass


def _selection():
    """
    A real RefGenomeSelection rather than a stub: resolve_input_genomes hands this
    straight to RunData.record_wanted_ref_tax_selection, which reads fields a stub
    silently wouldn't have.
    """
    attrition = PoolAttrition(resolved_rank="phylum", derep_rank=None, reps_only=True)
    return RefGenomeSelection(accessions=["GCA_000000001.1"], rows=[],
                              canonical="Testophyla", resolved_rank="phylum",
                              effective_derep_rank=None, warnings=[],
                              attrition=attrition)


@pytest.fixture
def captured_kwargs(monkeypatch):
    """Run resolve_input_genomes and hand back the kwargs the resolver was called with."""
    seen = {}

    def _fake_resolve(source, taxon, **kwargs):
        seen.update(kwargs)
        seen["source"] = source
        selection = _selection()
        return list(selection.accessions), selection

    import gtotree.utils.taxonomy.wanted_ref_tax as wrt
    monkeypatch.setattr(wrt, "resolve_wanted_ref_tax_accessions", _fake_resolve)
    monkeypatch.setattr(wrt, "expand_wanted_ref_tax",
                        lambda source, taxa: (list(taxa), []))
    # the source header reads the on-disk GTDB asset; not what's under test here
    monkeypatch.setattr("gtotree.utils.misc.general.wanted_ref_tax_source_line",
                        lambda *a, **kw: None)

    def _run(**arg_overrides):
        seen.clear()
        resolve_input_genomes(_Args(**arg_overrides), RunData(), _Error)
        return seen

    return _run


class TestPoolForwarding:

    def test_a_defaulted_gtdb_run_asks_for_representatives(self, captured_kwargs):
        assert captured_kwargs()["reps_only"] is True

    def test_gtdb_section_all_widens_the_pool(self, captured_kwargs):
        assert captured_kwargs(gtdb_section="all")["reps_only"] is False

    def test_gtdb_section_is_inert_under_source_ncbi(self, captured_kwargs):
        """
        None defers to NCBI's own default. Passing False here would look harmless and
        would quietly take ownership of a pool this flag doesn't govern.
        """
        seen = captured_kwargs(source="ncbi", gtdb_section="all")
        assert seen["reps_only"] is None


class TestMessagingVocabulary:

    def test_a_defaulted_run_does_not_blame_the_user_for_a_flag(self, captured_kwargs):
        """
        The pool is reps-only, but the user typed nothing. Reporting this as a
        user-specified filter sends them hunting for a flag they never passed.
        """
        seen = captured_kwargs()
        assert seen["reps_only"] is True
        assert seen["reps_only_requested"] is False

    def test_an_explicitly_typed_reps_may_be_blamed(self, captured_kwargs):
        assert captured_kwargs(gtdb_section="reps")["reps_only_requested"] is True

    def test_the_drivers_name_their_own_flag_not_dl_ncbi_assemblies(self, captured_kwargs):
        seen = captured_kwargs()
        assert seen["reps_flag"] == "--gtdb-section reps"
        assert "--representatives-only" not in str(seen["reps_flag"])

    def test_a_defaulted_run_still_says_how_to_widen(self, captured_kwargs):
        # the default branch used to dead-end at "this source's default pool" with no
        # way out; now there is one, so it has to be offered
        assert captured_kwargs()["reps_widen_hint"] == "--gtdb-section all"

    def test_gtdb_vocabulary_does_not_leak_into_an_ncbi_run(self, captured_kwargs):
        seen = captured_kwargs(source="ncbi")
        assert seen["reps_flag"] is None
        assert seen["reps_widen_hint"] is None
        assert seen["reps_only_requested"] is False
