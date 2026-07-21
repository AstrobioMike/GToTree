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
    d = dict(target_taxon=None, target_rank=None, derep_rank="off", source="refseq",
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
    _run(_args(target_taxon="Testophyla", source="refseq"))
    accs = _read_accs("ncbi-testophyla-phylum-refseq-accs.txt")
    assert sorted(accs) == ["GCF_000000001.1", "GCF_000000003.1"]  # only GCF


def test_source_scoping_precedes_derep(in_ncbi):
    """
    The regression this whole fix is about: --source refseq + derep must keep one GCF
    genome per family, NOT drop FamB because its best 'both'-pool genome was the GCA.
    """
    _run(_args(target_taxon="Testophyla", source="refseq", derep_rank="family"))
    accs = _read_accs("ncbi-testophyla-phylum-refseq-accs.txt")
    # one per family, and both must be GCF (the refseq-pool winners)
    assert len(accs) == 2
    assert all(a.startswith("GCF_") for a in accs)


def test_genbank_source_scopes_to_gca(in_ncbi):
    _run(_args(target_taxon="Testophyla", source="genbank", derep_rank="family"))
    accs = _read_accs("ncbi-testophyla-phylum-genbank-accs.txt")
    assert len(accs) == 2
    assert all(a.startswith("GCA_") for a in accs)


def test_both_source_derep_picks_best_regardless_of_prefix(in_ncbi):
    _run(_args(target_taxon="Testophyla", source="both", derep_rank="family"))
    accs = _read_accs("ncbi-testophyla-phylum-accs.txt")
    # FamB's best by quality is the 99.9-complete GCA_...004
    assert "GCA_000000004.1" in accs
    assert len(accs) == 2


def test_all_mode(in_ncbi):
    _run(_args(target_taxon="all", source="both"))
    accs = _read_accs("ncbi-all-accs.txt")
    assert len(accs) == 4


def test_taxid_mode(in_ncbi):
    _run(_args(target_taxon="5100", source="both"))  # FamA taxid
    accs = _read_accs("ncbi-taxid-5100-accs.txt")
    assert sorted(accs) == ["GCA_000000002.1", "GCF_000000001.1"]


def test_assembly_level_filter(in_ncbi):
    # only GCA_...004 is Scaffold; the rest Complete Genome
    _run(_args(target_taxon="Testophyla", source="both", assembly_level="scaffold"))
    accs = _read_accs("ncbi-testophyla-phylum-accs.txt")
    assert accs == ["GCA_000000004.1"]


def test_refseq_reference_genomes_only(in_ncbi):
    _run(_args(target_taxon="Testophyla", source="both",
               refseq_reference_genomes_only=True))
    accs = _read_accs("ncbi-testophyla-phylum-refseq-ref-accs.txt")
    # only GCF_...001 is a "reference genome"
    assert accs == ["GCF_000000001.1"]


def test_not_found_exits_cleanly(in_ncbi, capsys):
    code = _run(_args(target_taxon="Nonexistent", source="both"))
    assert code == 0
    assert "doesn't seem to exist" in capsys.readouterr().out


def test_bad_assembly_level_rejected(in_ncbi, capsys):
    code = _run(_args(target_taxon="Testophyla", assembly_level="banana"))
    assert code == 0
    assert "unrecognised" in capsys.readouterr().out.lower()


def test_coarser_derep_rank_rejected(in_ncbi, capsys):
    code = _run(_args(target_taxon="GenA", source="both", derep_rank="phylum"))
    assert code == 0
    assert capsys.readouterr().out  # friendly message emitted


def test_rank_counts(in_ncbi, capsys):
    _run(_args(get_rank_counts=True, source="refseq"))
    out = capsys.readouterr().out
    assert "phylum" in out
    assert "Unique Taxa" in out


def test_taxon_counts_is_case_insensitive(in_ncbi, capsys):
    # --get-taxon-counts routes through the shared resolver, so lowercase input
    # resolves to the canonical taxon rather than falling through to "no genomes"
    _run(_args(target_taxon="testophyla", source="both", get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "No genomes were found" not in out
    assert "Testophyla" in out                    # canonical casing echoed back


def test_taxon_counts_proper_case_no_match_note(in_ncbi, capsys):
    _run(_args(target_taxon="Testophyla", source="both", get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "Matched input" not in out


def test_taxon_counts_reports_per_rank_breakdown(in_ncbi, capsys):
    # GTDB-style format: "The rank 'X' has N <taxon> entries."
    _run(_args(target_taxon="Testophyla", source="both", get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "The rank 'phylum' has 4 Testophyla entries." in out


def test_taxon_counts_applies_source_filter(in_ncbi, capsys):
    # source refseq -> only the 2 GCF rows under Testophyla; scope note names the source
    _run(_args(target_taxon="Testophyla", source="refseq", get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "The rank 'phylum' has 2 Testophyla entries (in refseq)." in out


def test_taxon_counts_applies_assembly_level_filter(in_ncbi, capsys):
    # only GCA_...004 is Scaffold under Testophyla; scope note names the level
    _run(_args(target_taxon="Testophyla", source="both", assembly_level="scaffold",
               get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "The rank 'phylum' has 1 Testophyla entries (at assembly level Scaffold)." in out


def test_taxon_counts_reps_block(in_ncbi, capsys):
    # RefSeq-reference block, like GTDB's reps block
    _run(_args(target_taxon="Testophyla", source="both",
               refseq_reference_genomes_only=True, get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "The rank 'phylum' has 4 Testophyla entries." in out         # base pool (source both -> no scope note)
    assert "Of those, in considering only RefSeq reference genomes:" in out
    assert "has 1 Testophyla RefSeq reference genome entries." in out    # only GCF_...003


def test_taxon_counts_ignores_derep(in_ncbi, capsys):
    # --derep-rank must NOT collapse the count (derep is a pull-time reduction)
    _run(_args(target_taxon="Testophyla", source="both", derep_rank="family",
               get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "The rank 'phylum' has 4 Testophyla entries." in out   # all 4, not 2 families


def test_taxon_counts_multi_rank_breakdown(in_ncbi, capsys, tmp_path):
    # a name that appears at two ranks reported at both, no erroring
    records = [
        _rec("GCF_000000010.1", ("Bacteria", "Dualname", "ClassX", "OrdX", "FamX", "GenX", "GenX sp1")),
        _rec("GCF_000000011.1", ("Bacteria", "PhyY", "Dualname", "OrdY", "FamY", "GenY", "GenY sp1")),
    ]
    _write_mock_ncbi(tmp_path / PARQUET_FILENAME, records)
    _run(_args(target_taxon="Dualname", source="both", get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "The rank 'phylum' has 1 Dualname entries." in out
    assert "The rank 'class' has 1 Dualname entries." in out


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
