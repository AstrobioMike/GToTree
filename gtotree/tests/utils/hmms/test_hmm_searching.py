"""Unit tests for gtotree/main_stages/hmm_searching.py."""


from gtotree.utils.hmms.hmm_searching import (parse_hmmer_results,
                                                read_genome_hit_counts,
                                                rebuild_combined_SCG_outputs,
                                                write_genome_hit_counts)
from gtotree.utils.misc.general import RunData, SCGset
from gtotree.utils.misc.stages import GenomeRemovalStage


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
    """
    The dict is ordered by the target list, which is the order the per-genome counts
    file is written in. It's no longer what makes the counts *mean* anything -- that's
    the target name on each line now -- but a stable order keeps the file diffable.
    """
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
    from gtotree.utils.hmms import hmm_searching
    from gtotree.utils.misc.general import GenomeData

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
        from gtotree.utils.misc.general import RunData, GenomeData, SCGset

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
            gd.processing_done = True
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
        assert header.split("\t") == ["genome_id", "SCG_A", "SCG_B"]
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
        from gtotree.utils.misc.general import GenomeData
        rd = self._setup(tmp_path, ["G1", "G2"])
        ghost = GenomeData.from_acc("G_missing")
        ghost.processing_done = ghost.hmm_search_done = True
        rd.ncbi_accs.append(ghost)
        rd.update_all_input_genomes()
        rd.all_input_genomes[1].mark_removed("too few SCG hits",
                                         GenomeRemovalStage.SCG_HIT_FILTER)

        rebuild_combined_SCG_outputs(rd)
        _, rows = self._rows(rd)
        assert [r.split("\t")[0] for r in rows] == ["G1"]


# ---------------------------------------------------------------------------
# per-genome hit counts
# ---------------------------------------------------------------------------

class TestGenomeHitCountsFormat:
    """
    The per-genome counts file is keyed by target name.

    It used to be a bare positional row -- genome id, then one count per target in the
    order the target list had when *that genome* was searched -- while the combined
    table's header came from the target list as it stood at rebuild time. Those agree on
    a fresh run. They stop agreeing on a resume, because the NO_HITS check at the end of
    the search stage removes sets, so the second run's list is shorter: every count in
    every row written by the first run then sat under the wrong column name, and the row
    itself was wider than the header above it.
    """

    def _counts_file(self, tmp_path, text):
        path = tmp_path / "SCG-hit-counts.txt"
        path.write_text(text)
        return str(path)

    def test_a_round_trip_preserves_the_counts(self, tmp_path):
        path = str(tmp_path / "SCG-hit-counts.txt")
        write_genome_hit_counts(path, {"SCG_A": 1, "SCG_B": 0, "SCG_C": 3})

        assert read_genome_hit_counts(path, ["SCG_A", "SCG_B", "SCG_C"]) == {
            "SCG_A": 1, "SCG_B": 0, "SCG_C": 3}

    def test_counts_survive_the_target_list_shrinking_underneath_them(self, tmp_path):
        """
        The actual regression: written against three targets, read back after one was
        removed. The remaining two must still carry their own counts.
        """
        path = str(tmp_path / "SCG-hit-counts.txt")
        write_genome_hit_counts(path, {"SCG_A": 1, "SCG_B": 5, "SCG_C": 2})

        counts = read_genome_hit_counts(path, ["SCG_A", "SCG_C"])

        assert counts["SCG_A"] == 1
        assert counts["SCG_C"] == 2

    def test_a_legacy_row_is_read_when_its_width_still_matches(self, tmp_path):
        """
        A run already in flight when the format changed. Width matching is sufficient
        here: the target list only ever shrinks by removal and is never reordered, so an
        equal length means it is the same list.
        """
        path = self._counts_file(tmp_path, "G1\t1\t0\t3\n")

        assert read_genome_hit_counts(path, ["SCG_A", "SCG_B", "SCG_C"]) == {
            "SCG_A": 1, "SCG_B": 0, "SCG_C": 3}

    def test_a_legacy_row_of_the_wrong_width_is_refused_not_guessed(self, tmp_path):
        """
        This is the case that used to silently corrupt the table. Nothing in the old
        format records which targets those counts belonged to, so there is no way to
        recover the mapping -- None is the only honest answer.
        """
        path = self._counts_file(tmp_path, "G1\t1\t0\t3\n")

        assert read_genome_hit_counts(path, ["SCG_A", "SCG_C"]) is None

    def test_a_missing_file_reads_as_None(self, tmp_path):
        assert read_genome_hit_counts(str(tmp_path / "nope.txt"), ["SCG_A"]) is None

    def test_an_empty_file_reads_as_None(self, tmp_path):
        path = self._counts_file(tmp_path, "")
        assert read_genome_hit_counts(path, ["SCG_A"]) is None


class TestCombinedTableAlignment:

    def _setup(self, tmp_path, rows, scg_ids):
        """`rows` maps genome id -> the raw text of its counts file."""
        import os
        from gtotree.utils.misc.general import RunData, GenomeData, SCGset

        rd = RunData()
        rd.output_dir = str(tmp_path / "out")
        rd.hmm_results_dir = str(tmp_path / "hmm-results")
        rd.found_SCG_seqs_dir = str(tmp_path / "found-scgs")
        rd.general_ext = ".faa"
        for d in (rd.output_dir, rd.hmm_results_dir, rd.found_SCG_seqs_dir):
            os.makedirs(d, exist_ok=True)
        rd.SCG_targets = [SCGset.from_id(s) for s in scg_ids]

        for gid, text in rows.items():
            gd = GenomeData.from_acc(gid)
            gd.processing_done = True
            gd.hmm_search_done = True
            rd.ncbi_accs.append(gd)
            gdir = os.path.join(rd.hmm_results_dir, gid)
            os.makedirs(gdir, exist_ok=True)
            with open(os.path.join(gdir, "SCG-hit-counts.txt"), "w") as f:
                f.write(text)
        rd.update_all_input_genomes()
        return rd

    def _table(self, rd):
        import os
        path = os.path.join(rd.output_dir, "SCG-hit-counts.tsv")
        lines = [ln for ln in open(path).read().splitlines() if ln.strip()]
        header = lines[0].split("\t")
        return header, [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]

    def test_every_row_is_as_wide_as_the_header(self, tmp_path):
        """
        A row written against a longer target list used to be pasted in verbatim, so the
        table came out ragged -- unparseable by anything expecting a rectangle.
        """
        rd = self._setup(tmp_path,
                         {"G1": "target_SCG\tnum_hits\nSCG_A\t1\nSCG_B\t7\nSCG_C\t2\n"},
                         scg_ids=("SCG_A", "SCG_C"))

        rebuild_combined_SCG_outputs(rd)

        header, rows = self._table(rd)
        assert header == ["genome_id", "SCG_A", "SCG_C"]
        assert len(rows) == 1

    def test_counts_land_under_their_own_target(self, tmp_path):
        """
        SCG_B was removed between the two runs. Its count must not slide into SCG_C's
        column, which is exactly what the positional format did.
        """
        rd = self._setup(tmp_path,
                         {"G1": "target_SCG\tnum_hits\nSCG_A\t1\nSCG_B\t7\nSCG_C\t2\n"},
                         scg_ids=("SCG_A", "SCG_C"))

        rebuild_combined_SCG_outputs(rd)

        _header, rows = self._table(rd)
        assert rows[0] == {"genome_id": "G1", "SCG_A": "1", "SCG_C": "2"}

    def test_a_row_that_cannot_be_mapped_is_dropped_rather_than_misaligned(self, tmp_path):
        """
        A legacy row of the wrong width loses that genome from the counts table. Its
        sequences come from SCG-hits.faa and are unaffected, so the alignment and tree
        are unchanged -- a missing row is a far better outcome than a wrong one.
        """
        rd = self._setup(tmp_path,
                         {"G1": "G1\t1\t7\t2\n",
                          "G2": "target_SCG\tnum_hits\nSCG_A\t4\nSCG_C\t5\n"},
                         scg_ids=("SCG_A", "SCG_C"))

        rebuild_combined_SCG_outputs(rd)

        _header, rows = self._table(rd)
        assert [r["genome_id"] for r in rows] == ["G2"]

    def test_the_any_hit_tally_counts_the_right_targets(self, tmp_path):
        """
        The tally feeds `no_hits_reason`, so a misaligned read would attribute one set's
        hits to another and produce a removal reason naming the wrong gene.
        """
        rd = self._setup(tmp_path,
                         {"G1": "target_SCG\tnum_hits\nSCG_A\t0\nSCG_B\t7\nSCG_C\t2\n",
                          "G2": "target_SCG\tnum_hits\nSCG_A\t0\nSCG_B\t1\nSCG_C\t3\n"},
                         scg_ids=("SCG_A", "SCG_C"))

        rebuild_combined_SCG_outputs(rd)

        by_id = {s.id: s for s in rd.SCG_targets}
        assert by_id["SCG_A"].num_genomes_with_any_hit == 0
        assert by_id["SCG_C"].num_genomes_with_any_hit == 2
