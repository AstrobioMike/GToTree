"""Unit tests for gtotree/main_stages/hmm_searching.py."""


from gtotree.main_stages.hmm_searching import (parse_hmmer_results,
                                                rebuild_combined_SCG_outputs)
from gtotree.utils.general import RunData, SCGset


def _run_data(targets, best_hit_mode=False):
    rd = RunData()
    rd.SCG_targets = [SCGset.from_id(t) for t in targets]
    rd.best_hit_mode = best_hit_mode
    return rd


def _tblout(tmp_path, rows):
    """rows: (gene_id, target_SCG) in the order hmmsearch would have emitted them."""
    path = tmp_path / "SCG-hits-hmm.txt"
    with open(path, "w") as f:
        f.write("#  --- full sequence ----\n")
        f.write("# target name  accession  query name  accession  E-value\n")
        for gene, scg in rows:
            f.write(f"{gene}\t-\t{scg}\tACC\t1e-20\n")
    return str(path)


def test_single_hit_per_target(tmp_path):
    rd = _run_data(["A", "B"])
    path = _tblout(tmp_path, [("g1", "A"), ("g2", "B")])

    counts, gene_ids, n_hits, n_unique = parse_hmmer_results(path, rd)

    assert counts == {"A": 1, "B": 1}
    assert gene_ids == {"A": "g1", "B": "g2"}
    assert (n_hits, n_unique) == (2, 2)


def test_target_with_no_hits_counts_zero_and_has_no_gene(tmp_path):
    rd = _run_data(["A", "B", "C"])
    path = _tblout(tmp_path, [("g1", "A")])

    counts, gene_ids, n_hits, n_unique = parse_hmmer_results(path, rd)

    assert counts == {"A": 1, "B": 0, "C": 0}
    assert gene_ids["B"] is None and gene_ids["C"] is None
    assert (n_hits, n_unique) == (1, 1)


def test_multi_hit_counts_as_a_hit_but_not_unique(tmp_path):
    """Without --best-hit-mode a multi-hit SCG is counted but no sequence is pulled."""
    rd = _run_data(["A"], best_hit_mode=False)
    path = _tblout(tmp_path, [("g1", "A"), ("g2", "A"), ("g3", "A")])

    counts, gene_ids, n_hits, n_unique = parse_hmmer_results(path, rd)

    assert counts == {"A": 3}
    assert gene_ids["A"] is None
    assert (n_hits, n_unique) == (1, 0)


def test_best_hit_mode_takes_the_first_row(tmp_path):
    """hmmsearch emits best-scoring first, so the first row is the best hit."""
    rd = _run_data(["A"], best_hit_mode=True)
    path = _tblout(tmp_path, [("best", "A"), ("worse", "A"), ("worst", "A")])

    counts, gene_ids, n_hits, n_unique = parse_hmmer_results(path, rd)

    assert counts == {"A": 3}
    assert gene_ids["A"] == "best"
    assert (n_hits, n_unique) == (1, 0)


def test_best_hit_mode_picks_per_target_independently(tmp_path):
    """Interleaved targets must not let one SCG's first row leak into another's."""
    rd = _run_data(["A", "B"], best_hit_mode=True)
    path = _tblout(tmp_path, [
        ("a1", "A"), ("b1", "B"), ("a2", "A"), ("b2", "B"),
    ])

    _, gene_ids, _, _ = parse_hmmer_results(path, rd)

    assert gene_ids == {"A": "a1", "B": "b1"}


def test_hits_to_untargeted_scgs_are_ignored(tmp_path):
    """SCGs already dropped from the run contribute nothing to counts or totals."""
    rd = _run_data(["A"])
    path = _tblout(tmp_path, [("g1", "A"), ("g2", "REMOVED"), ("g3", "REMOVED")])

    counts, gene_ids, n_hits, n_unique = parse_hmmer_results(path, rd)

    assert counts == {"A": 1}
    assert set(gene_ids) == {"A"}
    assert (n_hits, n_unique) == (1, 1)


def test_counts_are_plain_ints(tmp_path):
    """Counts go straight into the per-genome counts file; numpy scalars would
    serialize differently."""
    rd = _run_data(["A"])
    path = _tblout(tmp_path, [("g1", "A")])
    counts, _, _, _ = parse_hmmer_results(path, rd)
    assert all(type(v) is int for v in counts.values())


def test_result_order_follows_the_target_list(tmp_path):
    """Counts are written positionally against the SCG-hit-counts.tsv header."""
    rd = _run_data(["A", "B", "C"])
    path = _tblout(tmp_path, [("g1", "C"), ("g2", "A")])
    counts, _, _, _ = parse_hmmer_results(path, rd)
    assert list(counts) == ["A", "B", "C"]


# ---------------------------------------------------------------------------
# worker contract
# ---------------------------------------------------------------------------

def test_worker_reports_failure_instead_of_raising(monkeypatch):
    """
    `run_pooled_stage` aborts a whole stage on a worker exception, so the worker reports
    problems in its status dict rather than propagating them.
    """
    from gtotree.main_stages import hmm_searching
    from gtotree.utils.general import GenomeData

    def exploding(*a, **k):
        raise MemoryError("blow-up inside the worker")

    monkeypatch.setattr(hmm_searching, "_hmm_search_worker_inner", exploding)

    status = hmm_searching._hmm_search_worker(GenomeData.from_acc("G1"), None)

    assert status["hmm_search_failed"] is True
    assert "MemoryError" in status.get("error", "")


# ---------------------------------------------------------------------------
# combined outputs
# ---------------------------------------------------------------------------

class TestRebuildCombinedSCGOutputs:
    """
    Builds SCG-hit-counts.tsv and the per-SCG FASTAs from each genome's own artifacts,
    so the result depends only on what is on disk, not on which genomes were searched in
    this invocation.
    """

    def _setup(self, tmp_path, genome_ids, scg_ids=("SCG_A", "SCG_B")):
        import os
        from gtotree.utils.general import RunData, GenomeData, SCGset

        rd = RunData()
        rd.output_dir = str(tmp_path / "out")
        rd.hmm_results_dir = str(tmp_path / "hmm-results")
        rd.found_SCG_seqs_dir = str(tmp_path / "found-scgs")
        rd.general_ext = ".faa"
        for d in (rd.output_dir, rd.hmm_results_dir, rd.found_SCG_seqs_dir):
            os.makedirs(d, exist_ok=True)
        rd.SCG_targets = [SCGset.from_id(s) for s in scg_ids]

        for gid in genome_ids:
            gd = GenomeData.from_acc(gid)
            gd.preprocessing_done = True
            gd.hmm_search_done = True
            rd.ncbi_accs.append(gd)
            gdir = os.path.join(rd.hmm_results_dir, gid)
            os.makedirs(gdir, exist_ok=True)
            with open(os.path.join(gdir, "SCG-hit-counts.txt"), "w") as f:
                f.write(f"{gid}\t1\t1\n")
            with open(os.path.join(gdir, "SCG-hits.faa"), "w") as f:
                for scg in scg_ids:
                    f.write(f">{scg}\nMSEQ{gid}\n")
        rd.update_all_input_genomes()
        return rd

    def _rows(self, rd):
        import os
        path = os.path.join(rd.output_dir, "SCG-hit-counts.tsv")
        lines = [ln for ln in open(path).read().splitlines() if ln.strip()]
        return lines[0], lines[1:]

    def test_includes_every_searched_genome(self, tmp_path):
        rd = self._setup(tmp_path, ["G1", "G2", "G3"])
        rebuild_combined_SCG_outputs(rd)
        header, rows = self._rows(rd)
        assert header.split("\t") == ["assembly_id", "SCG_A", "SCG_B"]
        assert [r.split("\t")[0] for r in rows] == ["G1", "G2", "G3"]

    def test_running_twice_produces_the_same_output(self, tmp_path):
        import os
        rd = self._setup(tmp_path, ["G1", "G2"])
        rebuild_combined_SCG_outputs(rd)
        first = (open(os.path.join(rd.output_dir, "SCG-hit-counts.tsv")).read(),
                 open(os.path.join(rd.found_SCG_seqs_dir, "SCG_A.faa")).read())
        rebuild_combined_SCG_outputs(rd)
        second = (open(os.path.join(rd.output_dir, "SCG-hit-counts.tsv")).read(),
                  open(os.path.join(rd.found_SCG_seqs_dir, "SCG_A.faa")).read())
        assert first == second
        assert first[1].count(">") == 2

    def test_overwrites_stale_content_rather_than_appending(self, tmp_path):
        import os
        rd = self._setup(tmp_path, ["G1", "G2"])
        with open(os.path.join(rd.found_SCG_seqs_dir, "SCG_A.faa"), "w") as f:
            f.write(">G1\nOLDJUNK\n")

        rebuild_combined_SCG_outputs(rd)
        content = open(os.path.join(rd.found_SCG_seqs_dir, "SCG_A.faa")).read()
        assert "OLDJUNK" not in content
        assert content.count(">G1") == 1

    def test_skips_genomes_with_no_artifacts_and_removed_genomes(self, tmp_path):
        from gtotree.utils.general import GenomeData
        rd = self._setup(tmp_path, ["G1", "G2"])
        ghost = GenomeData.from_acc("G_missing")
        ghost.preprocessing_done = ghost.hmm_search_done = True
        rd.ncbi_accs.append(ghost)
        rd.update_all_input_genomes()
        rd.all_input_genomes[1].mark_removed("too few SCG hits")

        rebuild_combined_SCG_outputs(rd)
        _, rows = self._rows(rd)
        assert [r.split("\t")[0] for r in rows] == ["G1"]
