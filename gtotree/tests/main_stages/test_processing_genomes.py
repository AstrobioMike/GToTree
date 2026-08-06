"""
Unit tests for gtotree/main_stages/processing_genomes.py.

The stage downloads/preprocesses and searches each genome in one worker call, then drops
that genome's sequence files. The tests below cover the observable results of that:
which genomes get processed, and what is left on disk afterwards.
"""

import datetime
import os
import shutil

import pytest  # type: ignore

from gtotree.tests.paths import DATA_DIR
from gtotree.utils.general import GenomeData
from gtotree.main_stages.processing_genomes import (SearchPlan,
                                                    genomes_needing_processing)

REQUIRED = ["prodigal"]
_missing = [b for b in REQUIRED if shutil.which(b) is None]


def _plan(do_pfam=False, do_ko=False, keep=False):
    return SearchPlan(do_pfam=do_pfam, do_ko=do_ko, keep_genome_files=keep)


def _genome(gid, **flags):
    gd = GenomeData.from_acc(gid)
    for k, v in flags.items():
        setattr(gd, k, v)
    return gd


class TestGenomesNeedingProcessing:
    """
    A genome's sequence files are deleted once its worker returns, so only fully
    processed genomes can be skipped -- a preprocessed-but-unsearched genome has nothing
    left on disk to search and is redone from the top.
    """

    def test_selects_unprocessed_and_partially_processed_genomes(self):
        done = _genome("G1", preprocessing_done=True, hmm_search_done=True)
        partial = _genome("G2", preprocessing_done=True, hmm_search_done=False)
        fresh = _genome("G3")
        removed = _genome("G4", removed=True)

        assert genomes_needing_processing([done, partial, fresh, removed], _plan()) == \
            [partial, fresh]

    def test_extra_searches_are_required_only_when_planned(self):
        gd = _genome("G1", preprocessing_done=True, hmm_search_done=True,
                     pfam_search_done=False)
        assert genomes_needing_processing([gd], _plan(do_pfam=True)) == [gd]
        assert genomes_needing_processing([gd], _plan(do_pfam=False)) == []


# ---------------------------------------------------------------------------
# end-to-end behaviour of the stage
# ---------------------------------------------------------------------------

def _process(tmp_path, source_flag, listing_paths, hmm, extra=()):
    from gtotree.cli.parser import parser
    from gtotree.utils.preflight_checks import preflight_checks
    from gtotree.main_stages.processing_genomes import process_genomes
    from gtotree.utils.context import log_file_var

    listing = tmp_path / "inputs.txt"
    listing.write_text("\n".join(str(p) for p in listing_paths) + "\n")

    argv = [source_flag, str(listing), "-H", str(hmm),
            "-o", str(tmp_path / "out"), "-j", "2", *extra]

    token = log_file_var.set(log_file_var.get())
    try:
        args, run_data = preflight_checks(parser().parse_args(argv))
        run_data.start_time = datetime.datetime.now()
        return args, process_genomes(args, run_data)
    finally:
        log_file_var.reset(token)


AA_GENOMES = [DATA_DIR / f"mock-{i}.faa" for i in range(1, 5)]
MOCK_HMM = DATA_DIR / "mock.hmm"


def test_processing_searches_every_genome_and_drops_its_sequence_files(tmp_path):
    _, run_data = _process(tmp_path, "-A", AA_GENOMES, MOCK_HMM)

    assert len(run_data.all_input_genomes) == 4
    assert all(gd.hmm_search_done for gd in run_data.all_input_genomes)
    assert os.listdir(run_data.ready_genome_files_dir) == []


def test_debug_keeps_the_sequence_files(tmp_path):
    _, run_data = _process(tmp_path, "-A", AA_GENOMES, MOCK_HMM, extra=["--debug"])
    assert sorted(os.listdir(run_data.ready_genome_files_dir)) == \
        ["mock-1.faa", "mock-2.faa", "mock-3.faa", "mock-4.faa"]


def test_running_the_stage_again_skips_everything_and_keeps_outputs_complete(tmp_path):
    from gtotree.main_stages.processing_genomes import process_genomes

    args, run_data = _process(tmp_path, "-A", AA_GENOMES, MOCK_HMM)
    table = os.path.join(run_data.output_dir, "SCG-hit-counts.tsv")
    first = open(table).read()

    process_genomes(args, run_data)

    assert open(table).read() == first
    assert len([ln for ln in first.splitlines() if ln.strip()]) == 5  # header + 4


# ---------------------------------------------------------------------------
# nucleotide mode
# ---------------------------------------------------------------------------

@pytest.mark.skipif(bool(_missing), reason=f"missing binaries: {_missing}")
def test_nucleotide_input_has_genes_called_and_targets_found(tmp_path, nt_genomes, nt_hmm):
    _, run_data = _process(tmp_path, "-f", nt_genomes, nt_hmm)

    assert all(gd.prodigal_used for gd in run_data.all_input_genomes)
    assert all(gd.num_genes > 50 for gd in run_data.all_input_genomes)
    assert all(gd.num_unique_SCG_hits == 1 for gd in run_data.all_input_genomes)
    # prodigal output must have its stop characters stripped
    for name in os.listdir(run_data.hmm_results_dir):
        hits = os.path.join(run_data.hmm_results_dir, name, "SCG-hits.faa")
        if os.path.isfile(hits):
            assert "*" not in open(hits).read()


@pytest.mark.skipif(bool(_missing), reason=f"missing binaries: {_missing}")
def test_nucleotide_mode_yields_nucleotide_hits_and_drops_both_files(
        tmp_path, nt_genomes, nt_hmm):
    """In `-z` a genome has both a .faa and a .fasta, and both are dropped."""
    _, run_data = _process(tmp_path, "-f", nt_genomes, nt_hmm, extra=["-z"])

    assert run_data.general_ext == ".fasta"
    hits = os.path.join(run_data.found_SCG_seqs_dir, "NT_MOCK.fasta")
    seqs, cur = {}, None
    for line in open(hits):
        if line.startswith(">"):
            cur = line[1:].strip(); seqs[cur] = ""
        else:
            seqs[cur] += line.strip()

    assert len(seqs) == 4
    for name, seq in seqs.items():
        assert set(seq) <= set("ACGTN"), f"{name} is not nucleotide"
        assert len(seq) % 3 == 0

    assert os.listdir(run_data.ready_genome_files_dir) == []
