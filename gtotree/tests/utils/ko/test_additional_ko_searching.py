"""
Unit tests for gtotree/utils/ko/additional_ko_searching.py.

The KO path is the one search in GToTree that shells out to an external binary
(kofamscan's `exec_annotation`) rather than going through pyhmmer, which is why it
needs this: the Pfam twin's failure handling is exercised incidentally by the
target-search integration tests, and the KO one was not exercised at all.

The failure branches are the point. A per-genome search that dies has to be recorded
and reported as failed, never raised -- a worker exception aborts the whole pooled
stage, taking every other genome's results down with it.

`exec_annotation` is never actually invoked; subprocess.run is patched throughout.
"""

import os
import subprocess
import types

import pytest  # type: ignore

from gtotree.utils.ko import additional_ko_searching as ko
from gtotree.utils.misc.general import SEARCH_FAILURE_FILENAME


MOCK_FAA = """\
>gene_1
MKAL
>gene_2
MTTV
>gene_3
MQQP
"""


def _failure_text(base_outpath):
    path = os.path.join(base_outpath, SEARCH_FAILURE_FILENAME)
    return open(path).read() if os.path.isfile(path) else None


# ---------------------------------------------------------------------------
# run_ko_search -- the external-binary call and its two failure modes
# ---------------------------------------------------------------------------

class TestRunKoSearch:

    def _run(self, tmp_path, monkeypatch, side_effect):
        base = str(tmp_path / "out") + "/"
        os.makedirs(base, exist_ok=True)
        monkeypatch.setattr(ko.subprocess, "run", side_effect)
        failed = ko.run_ko_search("profiles", "kos.tsv", base, "genome.faa")
        return failed, base

    def test_a_clean_run_reports_no_failure(self, tmp_path, monkeypatch):
        failed, base = self._run(
            tmp_path, monkeypatch,
            lambda cmd, **kw: subprocess.CompletedProcess(cmd, 0))

        assert failed is False
        assert _failure_text(base) is None

    def test_a_nonzero_exit_is_recorded_with_its_stderr(self, tmp_path, monkeypatch):
        def _fail(cmd, **kw):
            raise subprocess.CalledProcessError(3, cmd, stderr=b"profiles not found\n")

        failed, base = self._run(tmp_path, monkeypatch, _fail)

        assert failed is True
        text = _failure_text(base)
        assert "kofamscan failed" in text
        assert "exit 3" in text
        assert "profiles not found" in text

    def test_undecodable_stderr_does_not_take_the_search_down(self, tmp_path,
                                                              monkeypatch):
        # stderr is decoded with errors="replace"; invalid bytes must not raise
        def _fail(cmd, **kw):
            raise subprocess.CalledProcessError(1, cmd, stderr=b"\xff\xfe bad bytes")

        failed, base = self._run(tmp_path, monkeypatch, _fail)

        assert failed is True
        assert "exit 1" in _failure_text(base)

    def test_empty_stderr_is_handled(self, tmp_path, monkeypatch):
        def _fail(cmd, **kw):
            raise subprocess.CalledProcessError(1, cmd, stderr=None)

        failed, base = self._run(tmp_path, monkeypatch, _fail)

        assert failed is True
        assert "exit 1" in _failure_text(base)

    def test_a_missing_binary_is_recorded_rather_than_raised(self, tmp_path,
                                                             monkeypatch):
        """
        kofamscan is an optional dependency reached only via `-K`, so a user without it
        should get a recorded per-genome failure, not an OSError out of a worker.
        """
        def _boom(cmd, **kw):
            raise FileNotFoundError("no exec_annotation")

        failed, base = self._run(tmp_path, monkeypatch, _boom)

        assert failed is True
        assert "could not be run" in _failure_text(base)

    def test_the_tmp_dir_is_cleaned_up_on_success(self, tmp_path, monkeypatch):
        base = str(tmp_path / "out") + "/"
        os.makedirs(base + "kofamscan-tmp", exist_ok=True)
        monkeypatch.setattr(ko.subprocess, "run",
                            lambda cmd, **kw: subprocess.CompletedProcess(cmd, 0))

        ko.run_ko_search("profiles", "kos.tsv", base, "genome.faa")

        assert not os.path.exists(os.path.join(base, "kofamscan-tmp"))

    def test_the_tmp_dir_is_cleaned_up_on_failure(self, tmp_path, monkeypatch):
        # kofamscan's tmp dirs are large; leaving one per failed genome adds up
        base = str(tmp_path / "out") + "/"
        os.makedirs(base + "kofamscan-tmp", exist_ok=True)

        def _fail(cmd, **kw):
            raise subprocess.CalledProcessError(1, cmd, stderr=b"")
        monkeypatch.setattr(ko.subprocess, "run", _fail)

        ko.run_ko_search("profiles", "kos.tsv", base, "genome.faa")

        assert not os.path.exists(os.path.join(base, "kofamscan-tmp"))

    def test_the_command_carries_the_wanted_profiles_and_output(self, tmp_path,
                                                                monkeypatch):
        seen = {}

        def _capture(cmd, **kw):
            seen["cmd"] = cmd
            return subprocess.CompletedProcess(cmd, 0)
        monkeypatch.setattr(ko.subprocess, "run", _capture)

        base = str(tmp_path / "out") + "/"
        os.makedirs(base, exist_ok=True)
        ko.run_ko_search("/profiles", "/kos.tsv", base, "/genome.faa")

        cmd = seen["cmd"]
        assert cmd[0] == "exec_annotation"
        assert "/profiles" in cmd and "/kos.tsv" in cmd and "/genome.faa" in cmd
        # mapper format with unannotated rows dropped is what get_ko_counts parses
        assert "mapper" in cmd
        assert "--no-report-unannotated" in cmd


# ---------------------------------------------------------------------------
# _ko_search_worker -- exception containment
# ---------------------------------------------------------------------------

class TestWorkerContainment:

    def test_a_raising_worker_becomes_a_failure_dict(self, monkeypatch):
        """
        The wrapper exists because an exception escaping a pooled worker aborts the
        stage for every genome, not just this one.
        """
        def _boom(genome, run_data, aa_path=None):
            raise RuntimeError("something went sideways")
        monkeypatch.setattr(ko, "_ko_search_worker_inner", _boom)

        result = ko._ko_search_worker(types.SimpleNamespace(id="G1"), None)

        assert result["ko_search_failed"] is True
        assert "RuntimeError" in result["error"]
        assert "something went sideways" in result["error"]

    def test_a_clean_worker_result_passes_through(self, monkeypatch):
        monkeypatch.setattr(ko, "_ko_search_worker_inner",
                            lambda g, rd, aa_path=None: {"ko_search_failed": False})

        assert ko._ko_search_worker(types.SimpleNamespace(id="G1"), None) == \
            {"ko_search_failed": False}


# ---------------------------------------------------------------------------
# get_ko_counts
# ---------------------------------------------------------------------------

class TestGetKoCounts:

    def _results(self, tmp_path, body):
        path = tmp_path / "kofamscan-results.tsv"
        path.write_text(body)
        return str(path)

    def test_a_missing_results_file_is_all_zeros(self, tmp_path):
        # a genome whose search failed still needs a row in the counts table
        assert ko.get_ko_counts(["K00001", "K00002"],
                                str(tmp_path / "nope.tsv")) == [0, 0]

    def test_hits_are_counted_per_ko(self, tmp_path):
        path = self._results(tmp_path,
                             "gene_1\tK00001\ngene_2\tK00001\ngene_3\tK00002\n")
        assert ko.get_ko_counts(["K00001", "K00002"], path) == [2, 1]

    def test_the_order_follows_the_requested_kos(self, tmp_path):
        # the counts become columns in ko-hit-counts.tsv, so order has to be the
        # caller's, not the file's
        path = self._results(tmp_path, "gene_1\tK00002\ngene_2\tK00001\n")
        assert ko.get_ko_counts(["K00001", "K00002"], path) == [1, 1]
        assert ko.get_ko_counts(["K00002", "K00001"], path) == [1, 1]

    def test_a_ko_with_no_hits_is_zero_not_absent(self, tmp_path):
        path = self._results(tmp_path, "gene_1\tK00001\n")
        assert ko.get_ko_counts(["K00001", "K99999"], path) == [1, 0]

    def test_kos_outside_the_wanted_set_are_ignored(self, tmp_path):
        path = self._results(tmp_path, "gene_1\tK00001\ngene_2\tK55555\n")
        assert ko.get_ko_counts(["K00001"], path) == [1]

    def test_short_lines_are_skipped(self, tmp_path):
        # `-f mapper` emits a bare gene id for an unannotated gene
        path = self._results(tmp_path, "gene_1\ngene_2\tK00001\n\n")
        assert ko.get_ko_counts(["K00001"], path) == [1]

    def test_an_empty_results_file_is_all_zeros(self, tmp_path):
        assert ko.get_ko_counts(["K00001"], self._results(tmp_path, "")) == [0]


# ---------------------------------------------------------------------------
# the hit-sequence files
# ---------------------------------------------------------------------------

class TestWriteOutTmpKoHits:

    @pytest.fixture
    def faa(self, tmp_path):
        path = tmp_path / "genome.faa"
        path.write_text(MOCK_FAA)
        return str(path)

    def test_each_ko_gets_its_own_fasta(self, tmp_path, faa):
        results = tmp_path / "results.tsv"
        results.write_text("gene_1\tK00001\ngene_2\tK00002\n")
        out_base = str(tmp_path / "tmp") + "/"

        ko.write_out_tmp_ko_hits(str(results), out_base, faa)

        assert sorted(os.listdir(out_base)) == ["K00001.faa", "K00002.faa"]
        assert ">gene_1" in open(os.path.join(out_base, "K00001.faa")).read()

    def test_several_genes_hitting_one_ko_land_together(self, tmp_path, faa):
        results = tmp_path / "results.tsv"
        results.write_text("gene_1\tK00001\ngene_3\tK00001\n")
        out_base = str(tmp_path / "tmp") + "/"

        ko.write_out_tmp_ko_hits(str(results), out_base, faa)

        text = open(os.path.join(out_base, "K00001.faa")).read()
        assert text.count(">") == 2

    def test_a_gene_missing_from_the_fasta_is_skipped(self, tmp_path, faa):
        # rather than raising: the counts table is built from the results tsv, so a
        # desync here should cost a sequence, not the run
        results = tmp_path / "results.tsv"
        results.write_text("gene_1\tK00001\nnot_a_gene\tK00001\n")
        out_base = str(tmp_path / "tmp") + "/"

        ko.write_out_tmp_ko_hits(str(results), out_base, faa)

        assert open(os.path.join(out_base, "K00001.faa")).read().count(">") == 1


class TestCombineAllKoHits:

    def test_per_genome_files_are_concatenated_per_ko(self, tmp_path):
        tmp_area = tmp_path / "tmp"
        for genome, seq in (("G1", "MKAL"), ("G2", "MTTV")):
            d = tmp_area / genome
            d.mkdir(parents=True)
            (d / "K00001.faa").write_text(f">{genome}_gene\n{seq}\n")
        out_base = tmp_path / "out"
        out_base.mkdir()

        ko.combine_all_ko_hits(["K00001"], str(tmp_area), str(out_base))

        combined = (out_base / "K00001.faa").read_text()
        assert combined.count(">") == 2
        assert "G1_gene" in combined and "G2_gene" in combined

    def test_a_ko_nothing_hit_produces_no_file(self, tmp_path):
        tmp_area = tmp_path / "tmp"
        tmp_area.mkdir()
        out_base = tmp_path / "out"
        out_base.mkdir()

        ko.combine_all_ko_hits(["K00001"], str(tmp_area), str(out_base))

        assert os.listdir(out_base) == []


class TestWriteOutFailedKoTargets:

    def test_failed_targets_are_listed_one_per_line(self, tmp_path):
        run_data = types.SimpleNamespace(failed_ko_targets=["K00001", "K00002"],
                                         ko_results_dir=str(tmp_path))

        ko.write_out_failed_ko_targets(run_data)

        assert (tmp_path / "failed-ko-targets.txt").read_text() == "K00001\nK00002\n"

    def test_no_file_is_written_when_nothing_failed(self, tmp_path):
        run_data = types.SimpleNamespace(failed_ko_targets=[],
                                         ko_results_dir=str(tmp_path))

        ko.write_out_failed_ko_targets(run_data)

        assert not (tmp_path / "failed-ko-targets.txt").exists()


# ---------------------------------------------------------------------------
# write_ko_counts_table -- the matrix a user actually reads
# ---------------------------------------------------------------------------

class TestWriteKoCountsTable:

    def _run_data(self, tmp_path, genomes, targets):
        """genomes: list of (id, num_genes, ko_search_done, removed, {ko: n_hits})."""
        results_root = tmp_path / "individual-genome-results"
        gds = []
        for gid, num_genes, done, removed, hits in genomes:
            d = results_root / gid
            d.mkdir(parents=True)
            lines = []
            for ko_id, n in hits.items():
                lines += [f"{gid}_gene{i}\t{ko_id}" for i in range(n)]
            (d / "kofamscan-results.tsv").write_text("\n".join(lines) + "\n")
            gds.append(types.SimpleNamespace(id=gid, num_genes=num_genes,
                                             ko_search_done=done, removed=removed))
        return types.SimpleNamespace(all_input_genomes=gds,
                                     found_ko_targets=targets,
                                     ko_results_dir=str(tmp_path))

    def _table(self, tmp_path):
        return [ln.split("\t")
                for ln in (tmp_path / "ko-hit-counts.tsv").read_text().splitlines()]

    def test_the_header_is_the_targets_in_order(self, tmp_path):
        rd = self._run_data(tmp_path, [("G1", 10, True, False, {"K00001": 1})],
                            ["K00001", "K00002"])

        ko.write_ko_counts_table(rd)

        assert self._table(tmp_path)[0] == ["genome_id", "total_gene_count",
                                            "K00001", "K00002"]

    def test_counts_land_in_the_right_columns(self, tmp_path):
        rd = self._run_data(tmp_path,
                            [("G1", 10, True, False, {"K00001": 2, "K00002": 1})],
                            ["K00001", "K00002"])

        ko.write_ko_counts_table(rd)

        assert self._table(tmp_path)[1] == ["G1", "10", "2", "1"]

    def test_rows_are_sorted_by_genome_id(self, tmp_path):
        # a stable row order keeps the file diffable across runs
        rd = self._run_data(tmp_path, [("G2", 5, True, False, {"K00001": 1}),
                                       ("G1", 7, True, False, {"K00001": 1})],
                            ["K00001"])

        ko.write_ko_counts_table(rd)

        assert [r[0] for r in self._table(tmp_path)[1:]] == ["G1", "G2"]

    def test_removed_and_unsearched_genomes_are_left_out(self, tmp_path):
        rd = self._run_data(tmp_path, [("G1", 5, True, False, {"K00001": 1}),
                                       ("G2", 5, True, True, {"K00001": 1}),
                                       ("G3", 5, False, False, {"K00001": 1})],
                            ["K00001"])

        ko.write_ko_counts_table(rd)

        assert [r[0] for r in self._table(tmp_path)[1:]] == ["G1"]

    def test_the_tallies_come_back_per_genome(self, tmp_path):
        rd = self._run_data(tmp_path,
                            [("G1", 10, True, False, {"K00001": 3, "K00002": 0})],
                            ["K00001", "K00002"])

        tallies = ko.write_ko_counts_table(rd)

        # three hits, but only one of the two targets was hit
        assert tallies["G1"] == (3, 1)

    def test_no_genomes_still_writes_a_header(self, tmp_path):
        rd = self._run_data(tmp_path, [], ["K00001"])

        ko.write_ko_counts_table(rd)

        assert self._table(tmp_path)[0] == ["genome_id", "total_gene_count", "K00001"]
