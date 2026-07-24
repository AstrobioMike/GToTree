import argparse
import io
import os
import threading
import time

import pytest # type: ignore

import gtotree.utils.hmms.gen_scg_hmms_genomes as mod
from gtotree.utils.hmms.gen_scg_hmms import GenSCGHMMsError
from gtotree.utils.hmms.gen_scg_hmms_genomes import (
    MAX_DOWNLOAD_THREADS,
    MISSED_NO_PROTEINS,
    TargetGenomeError,
    fetch_amino_acids_pooled,
    genome_id_from_protein_name,
    read_accessions_file,
    relabel_and_append,
)


################################################################################
# accessions file
################################################################################

def test_read_accessions_file(tmp_path):
    path = tmp_path / "accs.txt"
    path.write_text("GCF_000000001.1\nGCA_000000002.1\n")
    assert read_accessions_file(str(path)) == ["GCF_000000001.1", "GCA_000000002.1"]


def test_read_accessions_file_skips_blanks_and_comments(tmp_path):
    path = tmp_path / "accs.txt"
    path.write_text("# a comment\n\nGCF_000000001.1\n   \nGCA_000000002.1\n")
    assert read_accessions_file(str(path)) == ["GCF_000000001.1", "GCA_000000002.1"]


def test_read_accessions_file_dedupes_preserving_order(tmp_path):
    path = tmp_path / "accs.txt"
    path.write_text("B\nA\nB\nC\nA\n")
    assert read_accessions_file(str(path)) == ["B", "A", "C"]


def test_read_accessions_file_missing(tmp_path):
    with pytest.raises(TargetGenomeError, match="can't be found"):
        read_accessions_file(str(tmp_path / "nope.txt"))


def test_read_accessions_file_empty(tmp_path):
    path = tmp_path / "accs.txt"
    path.write_text("\n# only comments\n")
    with pytest.raises(TargetGenomeError, match="no accessions"):
        read_accessions_file(str(path))


################################################################################
# header handling
################################################################################

def test_relabel_and_append_renames_and_strips_stops(tmp_path):
    src = tmp_path / "in.faa"
    src.write_text(">orig1 some description\nMKV*\n>orig2\nMA\nRT*\n")

    buf = io.StringIO()
    n = relabel_and_append("GCF_000091665.1", str(src), buf)

    assert n == 2
    assert buf.getvalue() == (">GCF_000091665.1_1\nMKV\n"
                              ">GCF_000091665.1_2\nMART\n")


def test_relabel_and_append_raises_when_no_proteins(tmp_path):
    src = tmp_path / "in.faa"
    src.write_text("")
    with pytest.raises(TargetGenomeError, match=MISSED_NO_PROTEINS):
        relabel_and_append("g1", str(src), io.StringIO())


@pytest.mark.parametrize("genome_id", [
    "GCF_000091665.1",
    "GCA_000007325.1",
    "my_genome_2",       # id that itself ends in _<number>
    "genome_1",
    "abc",
])
def test_genome_id_round_trips(genome_id):
    """
    Only the trailing `_<n>` may be stripped. Splitting on the FIRST underscore would
    turn GCF_000091665.1_1 into "GCF", silently merging every genome into one bucket.
    """
    assert genome_id_from_protein_name(f"{genome_id}_1") == genome_id
    assert genome_id_from_protein_name(f"{genome_id}_1234") == genome_id


################################################################################
# threaded fetching
################################################################################

def _pool_args(num_jobs):
    return argparse.Namespace(num_jobs=num_jobs)


def test_pooled_fetch_runs_concurrently(tmp_path, monkeypatch):
    live = [0]
    peak = [0]
    lock = threading.Lock()

    def fake_fetch(acc, entry, work_dir, nucleotide_fallback=True):
        with lock:
            live[0] += 1
            peak[0] = max(peak[0], live[0])
        time.sleep(0.02)
        path = os.path.join(work_dir, f"{acc}.faa")
        with open(path, "w") as f:
            f.write(">x\nMK\n")
        with lock:
            live[0] -= 1
        return path, False

    monkeypatch.setattr(mod, "fetch_amino_acids", fake_fetch)

    accs = [f"GCF_{i:09d}.1" for i in range(16)]
    got = []
    fetch_amino_acids_pooled(accs, {a: {} for a in accs}, str(tmp_path),
                             args=_pool_args(8),
                             on_result=lambda *a: got.append(a))

    assert len(got) == 16
    assert peak[0] > 1, "downloads did not run concurrently"
    assert peak[0] <= 8


def test_pooled_fetch_respects_thread_cap(tmp_path, monkeypatch):
    """NCBI concurrency is capped regardless of -j, to stay a good citizen."""
    live = [0]
    peak = [0]
    lock = threading.Lock()

    def fake_fetch(acc, entry, work_dir, nucleotide_fallback=True):
        with lock:
            live[0] += 1
            peak[0] = max(peak[0], live[0])
        time.sleep(0.02)
        with lock:
            live[0] -= 1
        return None, False

    monkeypatch.setattr(mod, "fetch_amino_acids", fake_fetch)

    accs = [f"G{i}" for i in range(60)]
    fetch_amino_acids_pooled(accs, {a: {} for a in accs}, str(tmp_path),
                             args=_pool_args(100), on_result=lambda *a: None)

    assert peak[0] <= MAX_DOWNLOAD_THREADS


def test_pooled_fetch_applies_results_on_one_thread(tmp_path, monkeypatch):
    """
    The apply callback must be single-threaded: the caller appends to one shared
    combined fasta from it, and concurrent appends would interleave records.
    """
    def fake_fetch(acc, entry, work_dir, nucleotide_fallback=True):
        time.sleep(0.01)
        return None, False

    monkeypatch.setattr(mod, "fetch_amino_acids", fake_fetch)

    seen_threads = set()

    def on_result(*args):
        seen_threads.add(threading.current_thread().name)

    accs = [f"G{i}" for i in range(20)]
    fetch_amino_acids_pooled(accs, {a: {} for a in accs}, str(tmp_path),
                             args=_pool_args(8), on_result=on_result)

    assert len(seen_threads) == 1


def test_pooled_fetch_isolates_failures(tmp_path, monkeypatch):
    """One bad accession must not abort the run; it's reported and the rest continue."""
    def fake_fetch(acc, entry, work_dir, nucleotide_fallback=True):
        if acc == "BAD_TARGET":
            raise TargetGenomeError("download failed")
        if acc == "BAD_UNEXPECTED":
            raise RuntimeError("kaboom")
        return os.path.join(work_dir, f"{acc}.faa"), False

    monkeypatch.setattr(mod, "fetch_amino_acids", fake_fetch)

    accs = ["OK_1", "BAD_TARGET", "OK_2", "BAD_UNEXPECTED", "OK_3"]
    results = {}
    fetch_amino_acids_pooled(accs, {a: {} for a in accs}, str(tmp_path),
                             args=_pool_args(4),
                             on_result=lambda acc, p, prod, err: results.__setitem__(acc, err))

    assert len(results) == 5
    assert results["OK_1"] is None
    assert results["BAD_TARGET"] == "download failed"
    assert "unexpected failure" in results["BAD_UNEXPECTED"]


def test_pooled_fetch_reports_prodigal_use(tmp_path, monkeypatch):
    def fake_fetch(acc, entry, work_dir, nucleotide_fallback=True):
        return None, acc.endswith("2")

    monkeypatch.setattr(mod, "fetch_amino_acids", fake_fetch)

    accs = ["G1", "G2"]
    flags = {}
    fetch_amino_acids_pooled(accs, {a: {} for a in accs}, str(tmp_path),
                             args=_pool_args(2),
                             on_result=lambda acc, p, prod, err: flags.__setitem__(acc, prod))

    assert flags == {"G1": False, "G2": True}


def test_pooled_fetch_handles_empty_list(tmp_path):
    # must not raise or spin up a pool for nothing
    fetch_amino_acids_pooled([], {}, str(tmp_path), args=_pool_args(4),
                             on_result=lambda *a: pytest.fail("should not be called"))


def test_pooled_fetch_accepts_num_jobs_without_args(tmp_path, monkeypatch):
    """Callers without an args namespace can pass num_jobs= directly."""
    monkeypatch.setattr(mod, "fetch_amino_acids",
                        lambda acc, entry, wd, nucleotide_fallback=True: (None, False))
    got = []
    fetch_amino_acids_pooled(["G1"], {"G1": {}}, str(tmp_path), num_jobs=2,
                             on_result=lambda *a: got.append(a))
    assert len(got) == 1


################################################################################
# prodigal
################################################################################

def test_run_prodigal_missing_binary(tmp_path, monkeypatch):
    def boom(*args, **kwargs):
        raise FileNotFoundError()

    monkeypatch.setattr(mod.subprocess, "run", boom)

    with pytest.raises(TargetGenomeError, match="prodigal doesn't seem to be available"):
        mod.run_prodigal(str(tmp_path / "in.fna"), str(tmp_path / "out.faa"))


def test_run_prodigal_cleans_up_on_failure(tmp_path, monkeypatch):
    """A failed call must not leave a `.part` behind."""
    out = tmp_path / "out.faa"

    class Result:
        returncode = 1
        stderr = b"boom"

    monkeypatch.setattr(mod.subprocess, "run", lambda *a, **k: Result())

    with pytest.raises(TargetGenomeError):
        mod.run_prodigal(str(tmp_path / "in.fna"), str(out))

    assert not out.exists()
    assert not (tmp_path / "out.faa.part").exists()


################################################################################
# fetch_amino_acids dispatch (download mocked; no network)
################################################################################

def _patch_downloader(monkeypatch, behavior):
    """
    Replace the shared atomic downloader. `behavior(url, dest)` either writes the file
    or raises, standing in for NCBI's response.
    """
    monkeypatch.setattr(mod, "_download_and_unzip", lambda: behavior)


def test_fetch_prefers_protein_file(tmp_path, monkeypatch):
    """The NCBI protein file is the submitter's own gene calls, so it wins."""
    urls = []

    def behavior(url, dest):
        urls.append(url)
        with open(dest, "w") as f:
            f.write(">p1\nMK\n")

    _patch_downloader(monkeypatch, behavior)

    entry = {"base_link": "https://example.org/GCF_000000001.1_ASM1"}
    path, used_prodigal = mod.fetch_amino_acids("GCF_000000001.1", entry, str(tmp_path))

    assert used_prodigal is False
    assert len(urls) == 1
    assert urls[0].endswith("_protein.faa.gz")


def test_fetch_falls_back_to_prodigal(tmp_path, monkeypatch):
    """
    GenBank assemblies frequently have no protein file, so a failure there means
    download the nucleotides and call genes instead.
    """
    urls = []

    def behavior(url, dest):
        urls.append(url)
        if url.endswith("_protein.faa.gz"):
            raise RuntimeError("404")
        with open(dest, "w") as f:
            f.write(">c1\nATGAAA\n")

    _patch_downloader(monkeypatch, behavior)
    monkeypatch.setattr(mod, "run_prodigal",
                        lambda nt, aa: open(aa, "w").write(">pred\nMK\n"))

    entry = {"base_link": "https://example.org/GCA_000000001.1_ASM1"}
    path, used_prodigal = mod.fetch_amino_acids("GCA_000000001.1", entry, str(tmp_path))

    assert used_prodigal is True
    assert [u.rsplit("_", 1)[-1] for u in urls] == ["protein.faa.gz", "genomic.fna.gz"]


def test_fetch_removes_nucleotide_temp(tmp_path, monkeypatch):
    def behavior(url, dest):
        if url.endswith("_protein.faa.gz"):
            raise RuntimeError("404")
        with open(dest, "w") as f:
            f.write(">c1\nATGAAA\n")

    _patch_downloader(monkeypatch, behavior)
    monkeypatch.setattr(mod, "run_prodigal",
                        lambda nt, aa: open(aa, "w").write(">pred\nMK\n"))

    mod.fetch_amino_acids("GCA_1.1", {"base_link": "https://x/GCA_1.1_ASM"},
                          str(tmp_path))

    leftovers = [f for f in os.listdir(tmp_path) if f.endswith("_genomic.fna")]
    assert leftovers == []


def test_fetch_without_fallback_raises(tmp_path, monkeypatch):
    def behavior(url, dest):
        raise RuntimeError("404")

    _patch_downloader(monkeypatch, behavior)

    with pytest.raises(TargetGenomeError, match="download failed"):
        mod.fetch_amino_acids("GCA_1.1", {"base_link": "https://x/GCA_1.1_ASM"},
                              str(tmp_path), nucleotide_fallback=False)


def test_fetch_raises_without_base_link(tmp_path):
    with pytest.raises(TargetGenomeError, match="no resolvable download location"):
        mod.fetch_amino_acids("GCA_1.1", {"base_link": ""}, str(tmp_path))


def test_fetch_raises_on_na_base_link(tmp_path):
    with pytest.raises(TargetGenomeError, match="no resolvable download location"):
        mod.fetch_amino_acids("GCA_1.1", {"base_link": "na"}, str(tmp_path))


def test_fetch_raises_when_nucleotide_download_fails(tmp_path, monkeypatch):
    def behavior(url, dest):
        raise RuntimeError("gone")

    _patch_downloader(monkeypatch, behavior)

    with pytest.raises(TargetGenomeError, match="download failed"):
        mod.fetch_amino_acids("GCA_1.1", {"base_link": "https://x/GCA_1.1_ASM"},
                              str(tmp_path))
