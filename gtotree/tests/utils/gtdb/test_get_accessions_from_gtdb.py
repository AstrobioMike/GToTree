import os
import glob
import types
import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pytest # type: ignore
from unittest.mock import patch

from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.gtdb.get_gtdb_data import PARQUET_FILENAME, VERSION_FILENAME
import gtotree.utils.gtdb.get_accessions_from_gtdb as mod

MODPATH = "gtotree.utils.gtdb.get_accessions_from_gtdb"

# columns the selection core / helper read for GTDB
_EXTRA_COLS = ["gtdb_representative", "ncbi_refseq_category",
               "checkm2_completeness", "checkm2_contamination",
               "genome_size", "contig_count"]


def _rec(gb_acc, lineage, rep="t", refseq="", comp="99.0", cont="0.5",
         size="4000000", contigs="1"):
    d = {"ncbi_genbank_assembly_accession": gb_acc,
         "gtdb_representative": rep, "ncbi_refseq_category": refseq,
         "checkm2_completeness": comp, "checkm2_contamination": cont,
         "genome_size": size, "contig_count": contigs}
    for i, r in enumerate(RANKS):
        d[r] = lineage[i]
    return d


def _write_mock_gtdb(path, records):
    keys = ["ncbi_genbank_assembly_accession"] + list(RANKS) + _EXTRA_COLS
    cols = {k: [rec[k] for rec in records] for k in keys}
    pq.write_table(pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}),
                   str(path))


# a small tree: phylum Testophyla with 2 classes, each with 2 genomes
_RECORDS = [
    _rec("GCA_000000001.1", ("Bacteria", "Testophyla", "ClassA", "OrdA", "FamA", "GenA", "GenA sp1")),
    _rec("GCA_000000002.1", ("Bacteria", "Testophyla", "ClassA", "OrdA", "FamA", "GenA", "GenA sp2"),
         comp="80.0"),
    _rec("GCA_000000003.1", ("Bacteria", "Testophyla", "ClassB", "OrdB", "FamB", "GenB", "GenB sp1"),
         refseq="reference genome"),
    _rec("GCA_000000004.1", ("Bacteria", "Testophyla", "ClassB", "OrdB", "FamB", "GenB", "GenB sp2"),
         rep="f"),
]


@pytest.fixture
def in_gtdb(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    _write_mock_gtdb(tmp_path / PARQUET_FILENAME, _RECORDS)
    (tmp_path / VERSION_FILENAME).write_text("r220\n2024-04-24\n")
    # get_gtdb_data() would try to download; short-circuit to "already present"
    monkeypatch.chdir(tmp_path)
    with patch(f"{MODPATH}.get_gtdb_data", return_value=str(tmp_path)):
        yield tmp_path


def _args(**kw):
    d = dict(wanted_ref_tax=None, target_rank=None, derep_rank="off",
             get_rank_counts=False, get_taxon_counts=False,
             gtdb_representatives_only=False, refseq_reference_genomes_only=False,
             get_table=False)
    d.update(kw)
    return types.SimpleNamespace(**d)


def _run(args):
    try:
        mod.get_accessions_from_gtdb(args)
        return None
    except SystemExit as e:
        return e.code


def _read_accs(pattern):
    matches = glob.glob(pattern)
    assert len(matches) == 1, f"expected 1 file for {pattern}, got {matches}"
    return [l.strip() for l in open(matches[0]) if l.strip()]


def test_plain_taxon_pull_derep_off(in_gtdb):
    _run(_args(wanted_ref_tax="Testophyla"))
    accs = _read_accs("gtdb-testophyla-phylum-accs.txt")
    assert len(accs) == 4  # all 4 genomes
    assert os.path.exists("gtdb-testophyla-phylum-metadata.tsv")


def test_derep_by_class_keeps_one_per_class(in_gtdb):
    _run(_args(wanted_ref_tax="Testophyla", derep_rank="class"))
    accs = _read_accs("gtdb-testophyla-phylum-accs.txt")
    assert len(accs) == 2  # one best per class (ClassA, ClassB)
    # ClassA's best by quality is the 99.0-complete GCA_000000001.1
    assert "GCA_000000001.1" in accs


def test_gtdb_representatives_only_filters(in_gtdb):
    _run(_args(wanted_ref_tax="Testophyla", gtdb_representatives_only=True))
    accs = _read_accs("gtdb-testophyla-phylum-gtdb-rep-accs.txt")
    # GCA_000000004.1 has gtdb_representative="f" -> excluded
    assert "GCA_000000004.1" not in accs
    assert len(accs) == 3


def test_all_taxon_bulk_dump(in_gtdb):
    _run(_args(wanted_ref_tax="all"))
    accs = _read_accs("gtdb-arc-and-bac-accs.txt")
    assert len(accs) == 4


def test_rank_counts(in_gtdb, capsys):
    _run(_args(get_rank_counts=True))
    out = capsys.readouterr().out
    assert "phylum" in out
    assert "Num. Unique Taxa" in out


def test_not_found_taxon_exits_cleanly(in_gtdb, capsys):
    code = _run(_args(wanted_ref_tax="Nonexistent"))
    assert code == 0
    assert "doesn't seem to exist" in capsys.readouterr().out


def test_coarser_derep_rank_is_rejected(in_gtdb, capsys):
    # target at genus, derep at phylum (coarser) -> ValueError translated
    code = _run(_args(wanted_ref_tax="GenA", derep_rank="phylum"))
    assert code == 0
    assert capsys.readouterr().out  # some friendly message emitted


def test_reps_flags_mutually_exclusive(in_gtdb, capsys):
    code = _run(_args(wanted_ref_tax="Testophyla",
                      gtdb_representatives_only=True,
                      refseq_reference_genomes_only=True))
    assert code == 1
    assert "Only one of" in capsys.readouterr().out


def test_taxon_counts_reports_the_dereplicated_size(in_gtdb, capsys):
    # --get-taxon-counts used to ignore --derep-rank entirely, so the number it quoted
    # was not the number a pull would return
    _run(_args(wanted_ref_tax="Testophyla", derep_rank="class", get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "The rank 'phylum' has 4 Testophyla entries." in out   # matches, unchanged
    assert "Dereplicated at 'class', that would be 2 genome(s)." in out


def test_taxon_counts_derep_size_matches_what_a_pull_returns(in_gtdb, capsys):
    _run(_args(wanted_ref_tax="Testophyla", derep_rank="class", get_taxon_counts=True))
    reported = capsys.readouterr().out

    _run(_args(wanted_ref_tax="Testophyla", derep_rank="class"))
    accs = _read_accs("gtdb-testophyla-phylum-accs.txt")
    assert f"that would be {len(accs)} genome(s)." in reported


def test_taxon_counts_without_derep_rank_says_nothing_about_derep(in_gtdb, capsys):
    _run(_args(wanted_ref_tax="Testophyla", get_taxon_counts=True))
    assert "Dereplicated at" not in capsys.readouterr().out


def test_rank_counts_are_scoped_to_a_wanted_taxon(in_gtdb, capsys):
    # -w used to be ignored here, printing whole-database counts
    _run(_args(wanted_ref_tax="Testophyla", get_rank_counts=True))
    out = capsys.readouterr().out
    assert "under 'Testophyla'" in out
    assert "domain" not in out          # coarser than the taxon's own rank
    assert "class      2" in out


def test_rank_counts_without_a_taxon_still_covers_the_database(in_gtdb, capsys):
    _run(_args(get_rank_counts=True))
    out = capsys.readouterr().out
    assert "domain" in out
    assert "under" not in out


def test_refseq_reference_only_counts_are_actually_filtered(in_gtdb, capsys):
    # the reps block was built against a literal "RefSeq" that never matched the
    # "refseq" value in play, so it silently reported UNFILTERED counts
    _run(_args(wanted_ref_tax="Testophyla", refseq_reference_genomes_only=True,
               get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "has 1 Testophyla refseq reference genome entries." in out


def test_all_with_refseq_reference_only_writes_only_reference_genomes(in_gtdb):
    # same stale literal: this used to write EVERY accession into a file named for
    # reference genomes
    _run(_args(wanted_ref_tax="all", refseq_reference_genomes_only=True))
    accs = _read_accs("gtdb-arc-and-bac-refseq-rep-accs.txt")
    assert accs == ["GCA_000000003.1"]


def test_all_with_derep_rank_dereplicates_within_each_domain(in_gtdb, capsys):
    """
    'all' has no rank of its own to group beneath, so it used to reject --derep-rank
    (and before that, silently ignore it). It now runs the per-domain selection and
    merges -- the fixture is one domain with 2 classes.
    """
    _run(_args(wanted_ref_tax="all", derep_rank="class"))
    accs = _read_accs("gtdb-arc-and-bac-accs.txt")
    assert len(accs) == 2
    assert "Dereplicating within each domain (Bacteria)" in capsys.readouterr().out
    assert os.path.exists("gtdb-arc-and-bac-metadata.tsv")


def test_all_with_derep_rank_spans_every_domain_present(tmp_path, monkeypatch, capsys):
    """The domain list comes from the asset, so nothing gets left behind."""
    records = _RECORDS + [
        _rec("GCA_000000005.1", ("Archaea", "ArcPhy", "ArcCls", "ArcOrd", "ArcFam",
                                 "ArcGen", "ArcGen sp1")),
    ]
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    _write_mock_gtdb(tmp_path / PARQUET_FILENAME, records)
    (tmp_path / VERSION_FILENAME).write_text("r220\n2024-04-24\n")
    monkeypatch.chdir(tmp_path)
    with patch(f"{MODPATH}.get_gtdb_data", return_value=str(tmp_path)):
        _run(_args(wanted_ref_tax="all", derep_rank="domain"))
    accs = _read_accs("gtdb-arc-and-bac-accs.txt")
    assert len(accs) == 2      # one genome per domain
    assert "Archaea, Bacteria" in capsys.readouterr().out


def test_all_taxon_counts_report_the_dereplicated_size(in_gtdb, capsys):
    _run(_args(wanted_ref_tax="all", derep_rank="class", get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "Dereplicated within each domain, that would be 2 genome(s)." in out


def test_derep_rank_off_with_all_still_works(in_gtdb):
    _run(_args(wanted_ref_tax="all", derep_rank="off"))
    assert len(_read_accs("gtdb-arc-and-bac-accs.txt")) == 4
