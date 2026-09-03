import glob
import types
import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pytest # type: ignore
from unittest.mock import patch

from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.ncbi.get_ncbi_assembly_data import PARQUET_FILENAME, DATE_FILENAME
import gtotree.utils.ncbi.get_accessions_from_ncbi as mod

MODPATH = "gtotree.utils.ncbi.get_accessions_from_ncbi"

_EXTRA_COLS = ["organism_name", "taxid", "asm_name", "assembly_level",
               "refseq_category", "checkm_completeness", "checkm_contamination",
               "genome_size", "genome_size_ungapped", "contig_count"]
_TAXID_COLS = [f"{r}_taxid" for r in RANKS]


def _rec(acc, lineage, level="Complete Genome", refseq="", comp="99.0", cont="0.5",
         taxid="100", lineage_taxids=None):
    d = {"assembly_accession": acc, "organism_name": "Testus " + lineage[-1],
         "taxid": taxid, "asm_name": "ASM1", "assembly_level": level,
         "refseq_category": refseq, "checkm_completeness": comp,
         "checkm_contamination": cont, "genome_size": "4000000",
         "genome_size_ungapped": "4000000", "contig_count": "1"}
    for i, r in enumerate(RANKS):
        d[r] = lineage[i]
    lineage_taxids = lineage_taxids or {}
    for r in RANKS:
        d[f"{r}_taxid"] = lineage_taxids.get(r, "0")
    return d


def _write_mock_ncbi(path, records):
    keys = ["assembly_accession"] + list(RANKS) + _EXTRA_COLS + _TAXID_COLS
    cols = {k: [rec.get(k, "") for rec in records] for k in keys}
    pq.write_table(pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}),
                   str(path))


# A tree under phylum Testophyla, two families each with a GCF (refseq) and GCA
# (genbank) genome. The GCA in FamB is higher quality, so a "both"-pool derep would
# pick it for FamB -- which is exactly what a naive post-filter would then wrongly drop
# under --source refseq.
_RECORDS = [
    _rec("GCF_000000001.1", ("Bacteria", "Testophyla", "ClassA", "OrdA", "FamA", "GenA", "GenA sp1"),
         refseq="reference genome", lineage_taxids={"phylum": "5000", "family": "5100"}),
    _rec("GCA_000000002.1", ("Bacteria", "Testophyla", "ClassA", "OrdA", "FamA", "GenA", "GenA sp2"),
         comp="80.0", lineage_taxids={"phylum": "5000", "family": "5100"}),
    _rec("GCF_000000003.1", ("Bacteria", "Testophyla", "ClassB", "OrdB", "FamB", "GenB", "GenB sp1"),
         comp="90.0", lineage_taxids={"phylum": "5000", "family": "5200"}),
    _rec("GCA_000000004.1", ("Bacteria", "Testophyla", "ClassB", "OrdB", "FamB", "GenB", "GenB sp2"),
         comp="99.9", level="Scaffold", lineage_taxids={"phylum": "5000", "family": "5200"}),
]


@pytest.fixture
def in_ncbi(tmp_path, monkeypatch):
    monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
    _write_mock_ncbi(tmp_path / PARQUET_FILENAME, _RECORDS)
    (tmp_path / DATE_FILENAME).write_text("2026,01,05\n")
    monkeypatch.chdir(tmp_path)
    with patch(f"{MODPATH}.get_ncbi_assembly_data", return_value=None):
        yield tmp_path


def _args(**kw):
    d = dict(wanted_ref_tax=None, target_rank=None, derep_rank="off", ncbi_section="refseq",
             assembly_level=None, refseq_reference_genomes_only=False,
             get_rank_counts=False, get_taxon_counts=False, get_table=False)
    d.update(kw)
    return types.SimpleNamespace(**d)


def _run(args):
    try:
        mod.get_accessions_from_ncbi(args)
        return None
    except SystemExit as e:
        return e.code


def _read_accs(pattern):
    matches = glob.glob(pattern)
    assert len(matches) == 1, f"expected 1 file for {pattern}, got {matches}"
    return [l.strip() for l in open(matches[0]) if l.strip()]


def test_taxon_refseq_derep_off(in_ncbi):
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="refseq"))
    accs = _read_accs("ncbi-testophyla-phylum-refseq-accs.txt")
    assert sorted(accs) == ["GCF_000000001.1", "GCF_000000003.1"]  # only GCF


def test_source_scoping_precedes_derep(in_ncbi):
    """
    Source scoping applies before dereplication: with --source refseq, each family's
    representative is picked from the RefSeq pool only. Otherwise a family whose overall
    best genome is a GCA would be dropped entirely rather than represented by its best
    GCF.
    """
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="refseq", derep_rank="family"))
    accs = _read_accs("ncbi-testophyla-phylum-refseq-accs.txt")
    # one per family, and both must be GCF (the refseq-pool winners)
    assert len(accs) == 2
    assert all(a.startswith("GCF_") for a in accs)


def test_genbank_source_scopes_to_gca(in_ncbi):
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="genbank", derep_rank="family"))
    accs = _read_accs("ncbi-testophyla-phylum-genbank-accs.txt")
    assert len(accs) == 2
    assert all(a.startswith("GCA_") for a in accs)


def test_both_source_derep_picks_best_regardless_of_prefix(in_ncbi):
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both", derep_rank="family"))
    accs = _read_accs("ncbi-testophyla-phylum-accs.txt")
    # FamB's best by quality is the 99.9-complete GCA_...004
    assert "GCA_000000004.1" in accs
    assert len(accs) == 2


def test_all_mode(in_ncbi):
    _run(_args(wanted_ref_tax="all", ncbi_section="both"))
    accs = _read_accs("ncbi-all-accs.txt")
    assert len(accs) == 4


def test_taxid_mode(in_ncbi):
    _run(_args(wanted_ref_tax="5100", ncbi_section="both"))  # FamA taxid
    accs = _read_accs("ncbi-taxid-5100-accs.txt")
    assert sorted(accs) == ["GCA_000000002.1", "GCF_000000001.1"]


def test_assembly_level_filter(in_ncbi):
    # only GCA_...004 is Scaffold; the rest Complete Genome
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both", assembly_level="scaffold"))
    accs = _read_accs("ncbi-testophyla-phylum-accs.txt")
    assert accs == ["GCA_000000004.1"]


def test_refseq_reference_genomes_only(in_ncbi):
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both",
               refseq_reference_genomes_only=True))
    accs = _read_accs("ncbi-testophyla-phylum-refseq-ref-accs.txt")
    # only GCF_...001 is a "reference genome"
    assert accs == ["GCF_000000001.1"]


def test_not_found_exits_cleanly(in_ncbi, capsys):
    code = _run(_args(wanted_ref_tax="Nonexistent", ncbi_section="both"))
    assert code == 0
    assert "doesn't seem to exist" in capsys.readouterr().out


def test_bad_assembly_level_rejected(in_ncbi, capsys):
    code = _run(_args(wanted_ref_tax="Testophyla", assembly_level="banana"))
    assert code == 0
    assert "unrecognised" in capsys.readouterr().out.lower()


def test_cli_assembly_level_is_repeatable():
    # repeated --assembly-level must accumulate (argparse append), not overwrite down to
    # the last value
    from gtotree.utils.taxonomy.get_accs_shared import parse_assembly_levels
    args = mod.build_parser().parse_args(
        ["-w", "Testophyla", "--assembly-level", "complete",
         "--assembly-level", "contig"])
    assert args.assembly_level == ["complete", "contig"]
    assert parse_assembly_levels(args.assembly_level) == ["Complete Genome", "Contig"]


def test_coarser_derep_rank_rejected(in_ncbi, capsys):
    code = _run(_args(wanted_ref_tax="GenA", ncbi_section="both", derep_rank="phylum"))
    assert code == 0
    assert capsys.readouterr().out  # friendly message emitted


def test_rank_counts(in_ncbi, capsys):
    _run(_args(get_rank_counts=True, ncbi_section="refseq"))
    out = capsys.readouterr().out
    assert "phylum" in out
    assert "Unique Taxa" in out


def test_taxon_counts_is_case_insensitive(in_ncbi, capsys):
    # --get-taxon-counts routes through the shared resolver, so lowercase input
    # resolves to the canonical taxon rather than falling through to "no genomes"
    _run(_args(wanted_ref_tax="testophyla", ncbi_section="both", get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "No genomes were found" not in out
    assert "Testophyla" in out                    # canonical casing echoed back


def test_taxon_counts_proper_case_no_match_note(in_ncbi, capsys):
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both", get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "Matched input" not in out


def test_taxon_counts_reports_per_rank_breakdown(in_ncbi, capsys):
    # GTDB-style format, now with the generic filter tag appended to every count line
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both", get_taxon_counts=True))
    out = capsys.readouterr().out
    assert ("The rank 'phylum' has 4 Testophyla entries "
            "(after any specified filters)." in out)


def test_taxon_counts_applies_source_filter(in_ncbi, capsys):
    # source refseq -> only the 2 GCF rows under Testophyla; the count reflects the
    # filter, and the generic "(after any specified filters)" tag flags that it does
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="refseq", get_taxon_counts=True))
    out = capsys.readouterr().out
    assert ("The rank 'phylum' has 2 Testophyla entries "
            "(after any specified filters)." in out)


def test_taxon_counts_applies_assembly_level_filter(in_ncbi, capsys):
    # only GCA_...004 is Scaffold under Testophyla; count drops to 1 and carries the tag
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both", assembly_level="scaffold",
               get_taxon_counts=True))
    out = capsys.readouterr().out
    assert ("The rank 'phylum' has 1 Testophyla entries "
            "(after any specified filters)." in out)


def test_taxon_counts_reps_block(in_ncbi, capsys):
    # RefSeq-reference block, like GTDB's reps block
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both",
               refseq_reference_genomes_only=True, get_taxon_counts=True))
    out = capsys.readouterr().out
    assert ("The rank 'phylum' has 4 Testophyla entries "
            "(after any specified filters)." in out)  # base pool
    assert "Of those, in considering only RefSeq reference genomes:" in out
    assert "has 1 Testophyla RefSeq reference genome entries." in out    # only GCF_...003


def test_taxon_counts_reports_matches_and_the_dereplicated_size(in_ncbi, capsys):
    # the match count must NOT be collapsed by --derep-rank, but the dereplicated size
    # is reported alongside it -- it used to be omitted entirely
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both", derep_rank="family",
               get_taxon_counts=True))
    out = capsys.readouterr().out
    assert ("The rank 'phylum' has 4 Testophyla entries "
            "(after any specified filters)." in out)  # all 4, not 2 families
    assert "Dereplicated at 'family', that would be 2 genome(s)." in out


def test_taxon_counts_multi_rank_breakdown(in_ncbi, capsys, tmp_path):
    # a name that appears at two ranks reported at both, no erroring
    records = [
        _rec("GCF_000000010.1", ("Bacteria", "Dualname", "ClassX", "OrdX", "FamX", "GenX", "GenX sp1")),
        _rec("GCF_000000011.1", ("Bacteria", "PhyY", "Dualname", "OrdY", "FamY", "GenY", "GenY sp1")),
    ]
    _write_mock_ncbi(tmp_path / PARQUET_FILENAME, records)
    _run(_args(wanted_ref_tax="Dualname", ncbi_section="both", get_taxon_counts=True))
    out = capsys.readouterr().out
    assert ("The rank 'phylum' has 1 Dualname entries "
            "(after any specified filters)." in out)
    assert ("The rank 'class' has 1 Dualname entries "
            "(after any specified filters)." in out)


def test_get_table_writes_full_metadata_tsv(in_ncbi, capsys):
    import csv
    _run(_args(get_table=True))
    out = capsys.readouterr().out
    assert "NCBI table written to" in out
    matches = glob.glob("ncbi-assembly-summary-metadata.tsv")
    assert len(matches) == 1
    with open(matches[0]) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    # full dump includes columns beyond the selection subset (_COLUMNS); phylum_taxid
    # is written by the mock but is NOT in the selection subset, so it proves the dump
    # is the full table
    assert len(rows) == len(_RECORDS)
    assert "phylum_taxid" in rows[0]


def test_assembly_level_is_applied_before_dereplication(in_ncbi):
    """
    FamB holds a Complete Genome (GCF_...003) and a higher-quality Scaffold
    (GCA_...004). Filtering on level AFTER dereplication picked the Scaffold as FamB's
    best and then deleted it, so FamB contributed nothing at all -- instead of
    contributing the best genome that actually meets the requested level.
    """
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both", derep_rank="family",
               assembly_level=["complete"]))
    accs = _read_accs("ncbi-testophyla-phylum-accs.txt")
    assert sorted(accs) == ["GCF_000000001.1", "GCF_000000003.1"]


def test_assembly_level_still_filters_with_derep_off(in_ncbi):
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both", assembly_level=["scaffold"]))
    assert _read_accs("ncbi-testophyla-phylum-accs.txt") == ["GCA_000000004.1"]


def test_assembly_level_still_filters_for_all(in_ncbi):
    _run(_args(wanted_ref_tax="all", ncbi_section="both", assembly_level=["scaffold"]))
    assert _read_accs("ncbi-all-accs.txt") == ["GCA_000000004.1"]


def test_assembly_level_still_filters_for_a_taxid(in_ncbi):
    _run(_args(wanted_ref_tax="5000", ncbi_section="both", assembly_level=["scaffold"]))
    assert _read_accs("ncbi-taxid-5000-accs.txt") == ["GCA_000000004.1"]


def test_taxon_counts_derep_size_matches_what_a_pull_returns(in_ncbi, capsys):
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both", derep_rank="family",
               get_taxon_counts=True))
    reported = capsys.readouterr().out

    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both", derep_rank="family"))
    accs = _read_accs("ncbi-testophyla-phylum-accs.txt")
    assert f"that would be {len(accs)} genome(s)." in reported


def test_taxon_counts_derep_size_respects_the_source_filter(in_ncbi, capsys):
    # only 2 GCF rows under Testophyla, one per family -> 2 either way, but the
    # dereplicated number must come from the refseq-scoped pool
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="refseq", derep_rank="genus",
               get_taxon_counts=True))
    reported = capsys.readouterr().out

    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="refseq", derep_rank="genus"))
    accs = _read_accs("ncbi-testophyla-phylum-refseq-accs.txt")
    assert f"that would be {len(accs)} genome(s)." in reported


def test_taxon_counts_without_derep_rank_says_nothing_about_derep(in_ncbi, capsys):
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both", get_taxon_counts=True))
    assert "Dereplicated at" not in capsys.readouterr().out


def test_rank_counts_are_scoped_to_a_wanted_taxon(in_ncbi, capsys):
    # -w used to be ignored here, printing whole-database counts
    _run(_args(wanted_ref_tax="Testophyla", ncbi_section="both", get_rank_counts=True))
    out = capsys.readouterr().out
    assert "under 'Testophyla'" in out
    assert "domain" not in out          # coarser than the taxon's own rank
    assert "family     2" in out


def test_rank_counts_without_a_taxon_still_covers_the_database(in_ncbi, capsys):
    _run(_args(get_rank_counts=True, ncbi_section="both"))
    out = capsys.readouterr().out
    assert "domain" in out
    assert "under" not in out


def test_all_with_derep_rank_dereplicates_within_each_domain(in_ncbi, capsys):
    _run(_args(wanted_ref_tax="all", ncbi_section="both", derep_rank="family"))
    accs = _read_accs("ncbi-all-accs.txt")
    assert len(accs) == 2      # one per family, within the one domain present
    assert "Dereplicating within each domain (Bacteria)" in capsys.readouterr().out


def test_all_with_derep_rank_includes_eukaryotes(in_ncbi, tmp_path, capsys):
    """
    The domain list is read from the table, not hardcoded to the prokaryotic pair --
    NCBI carries eukaryotes and they must not be silently dropped.
    """
    records = _RECORDS + [
        _rec("GCF_000000020.1", ("Eukaryota", "Ascomycota", "Saccharomycetes",
                                 "Saccharomycetales", "Saccharomycetaceae",
                                 "Saccharomyces", "Saccharomyces cerevisiae")),
    ]
    _write_mock_ncbi(tmp_path / PARQUET_FILENAME, records)
    _run(_args(wanted_ref_tax="all", ncbi_section="both", derep_rank="domain"))
    accs = _read_accs("ncbi-all-accs.txt")
    assert len(accs) == 2      # one bacterium, one eukaryote
    assert "GCF_000000020.1" in accs
    assert "Bacteria, Eukaryota" in capsys.readouterr().out


def test_derep_rank_with_a_taxid_is_refused(in_ncbi, capsys):
    code = _run(_args(wanted_ref_tax="5000", ncbi_section="both", derep_rank="family"))
    assert code == 0
    assert "`--derep-rank` can't be applied with a taxid" in capsys.readouterr().out
    assert not glob.glob("ncbi-*accs.txt")


def test_derep_rank_off_with_a_taxid_still_works(in_ncbi):
    _run(_args(wanted_ref_tax="5000", ncbi_section="both"))
    assert len(_read_accs("ncbi-taxid-5000-accs.txt")) == 4


################################################################################
# 'all' and the genomes it can't reach
#
# The NCBI table carries rows with NO assigned domain (viral and
# metagenome/uncultured entries). 'all' is expanded to the table's DOMAINS, so those
# rows are unreachable by it -- which is correct for GToTree, but has to be both
# CONSISTENT (the same pool with and without --derep-rank) and REPORTED (an unscoped
# --get-rank-counts otherwise quotes a number no pull will produce).
################################################################################

# two viral rows: no domain, and a class that exists nowhere else in the table
_VIRAL_RECORDS = [
    _rec("GCF_000000030.1", ("NA", "Uroviricota", "Caudoviricetes", "NA", "NA", "NA",
                             "Phage sp1")),
    _rec("GCF_000000031.1", ("NA", "Uroviricota", "Caudoviricetes", "NA", "NA", "NA",
                             "Phage sp2")),
]


@pytest.fixture
def in_ncbi_with_viruses(tmp_path, monkeypatch):
    monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
    _write_mock_ncbi(tmp_path / PARQUET_FILENAME, _RECORDS + _VIRAL_RECORDS)
    (tmp_path / DATE_FILENAME).write_text("2026,01,05\n")
    monkeypatch.chdir(tmp_path)
    with patch(f"{MODPATH}.get_ncbi_assembly_data", return_value=None):
        yield tmp_path


def test_all_without_derep_excludes_domainless_genomes(in_ncbi_with_viruses):
    """
    'all' has to mean the same pool whether or not --derep-rank is set. It used to be
    a whole-table dump without derep (viruses in) and a domain walk with it (viruses
    out), so an unrelated flag decided whether viral genomes appeared.
    """
    _run(_args(wanted_ref_tax="all", ncbi_section="both"))
    accs = _read_accs("ncbi-all-accs.txt")
    assert len(accs) == 4
    assert not any(a in accs for a in ("GCF_000000030.1", "GCF_000000031.1"))


def test_all_reports_the_genomes_it_left_behind(in_ncbi_with_viruses, capsys):
    _run(_args(wanted_ref_tax="all", ncbi_section="both", derep_rank="class"))
    out = capsys.readouterr().out
    assert "2 genome(s)" in out
    assert "no assigned domain" in out


def test_rank_counts_with_all_are_scoped_to_what_all_pulls(in_ncbi_with_viruses, capsys):
    """
    The reconciliation: the class count reported for `-w all` must equal the number of
    accessions `-w all --derep-rank class` writes. Unscoped it counted the viral class
    too, which is the shape of "--get-rank-counts says 326 but we get 280".
    """
    _run(_args(wanted_ref_tax="all", ncbi_section="both", get_rank_counts=True))
    counts_out = capsys.readouterr().out
    assert "class      2" in counts_out          # ClassA + ClassB, not Caudoviricetes

    _run(_args(wanted_ref_tax="all", ncbi_section="both", derep_rank="class"))
    assert len(_read_accs("ncbi-all-accs.txt")) == 2


def test_rank_counts_without_a_taxon_still_counts_the_whole_table(in_ncbi_with_viruses,
                                                                 capsys):
    """Unscoped `--get-rank-counts` is documented as the whole database -- unchanged."""
    _run(_args(get_rank_counts=True, ncbi_section="both"))
    assert "class      3" in capsys.readouterr().out   # includes Caudoviricetes


def test_domainless_genomes_are_still_reachable_by_name(in_ncbi_with_viruses):
    """They're excluded from 'all', not from the tool."""
    _run(_args(wanted_ref_tax="Uroviricota", ncbi_section="both"))
    accs = _read_accs("ncbi-uroviricota-phylum-accs.txt")
    assert sorted(accs) == ["GCF_000000030.1", "GCF_000000031.1"]


def test_eukaryote_alias_resolves(in_ncbi, tmp_path):
    """`-w eukarya` (and friends) resolve to Eukaryota wherever a taxon is taken."""
    records = _RECORDS + [
        _rec("GCF_000000020.1", ("Eukaryota", "Ascomycota", "Saccharomycetes",
                                 "Saccharomycetales", "Saccharomycetaceae",
                                 "Saccharomyces", "Saccharomyces cerevisiae")),
    ]
    _write_mock_ncbi(tmp_path / PARQUET_FILENAME, records)
    _run(_args(wanted_ref_tax="eukarya", ncbi_section="both"))
    accs = _read_accs("ncbi-eukaryota-domain-accs.txt")
    assert accs == ["GCF_000000020.1"]
