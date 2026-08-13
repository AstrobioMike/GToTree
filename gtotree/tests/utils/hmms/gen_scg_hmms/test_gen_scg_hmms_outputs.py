import pytest # type: ignore
from gtotree.utils.hmms.gen_scg_hmms import gen_scg_hmms_outputs as outputs
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import PfamProfileInfo


def _read(path):
    with open(path) as f:
        return [line.rstrip("\n").split("\t") for line in f]


################################################################################
# SCG targets info
################################################################################

def test_write_scg_targets_info(tmp_path):
    info = {"PF00001.27": PfamProfileInfo("PF00001.27", "7tm_1", "A receptor", 66.58)}
    path = outputs.write_scg_targets_info(str(tmp_path), ["PF00001.27"], info)

    rows = _read(path)
    assert rows[0] == ["pfam_id", "name", "description", "average_coverage"]
    assert rows[1] == ["PF00001.27", "7tm_1", "A receptor", "66.58"]


def test_write_scg_targets_info_falls_back_to_na(tmp_path):
    """A retained profile with no info row must still appear, not vanish."""
    path = outputs.write_scg_targets_info(str(tmp_path), ["PF99999.1"], {})
    rows = _read(path)
    assert rows[1] == ["PF99999.1", "NA", "NA", "NA"]


def test_write_scg_targets_info_preserves_order(tmp_path):
    info = {a: PfamProfileInfo(a, a, "d", 60.0) for a in ["PF3.1", "PF1.1", "PF2.1"]}
    path = outputs.write_scg_targets_info(str(tmp_path), ["PF3.1", "PF1.1", "PF2.1"], info)
    assert [r[0] for r in _read(path)[1:]] == ["PF3.1", "PF1.1", "PF2.1"]


################################################################################
# hit counts
################################################################################

def test_write_hit_counts_matrix(tmp_path):
    path = outputs.write_hit_counts(
        str(tmp_path), ["G1", "G2"], ["PF1.1", "PF2.1"],
        {"G1": {"PF1.1": 1}, "G2": {"PF2.1": 2}})

    rows = _read(path)
    assert rows[0] == ["genome", "PF1.1", "PF2.1"]
    assert rows[1] == ["G1", "1", "0"]
    assert rows[2] == ["G2", "0", "2"]


def test_hit_counts_include_all_searched_profiles(tmp_path):
    """
    Columns cover every profile SEARCHED, not just the retained ones, so the table can
    show why a marker was dropped.
    """
    path = outputs.write_hit_counts(
        str(tmp_path), ["G1"], ["kept.1", "dropped.1"], {"G1": {"kept.1": 1}})
    assert _read(path)[0] == ["genome", "kept.1", "dropped.1"]


def test_hit_counts_row_per_genome_even_when_no_hits(tmp_path):
    path = outputs.write_hit_counts(str(tmp_path), ["G1", "G2"], ["PF1.1"],
                                    {"G1": {"PF1.1": 1}})
    rows = _read(path)
    assert len(rows) == 3          # header + 2 genomes
    assert rows[2] == ["G2", "0"]


################################################################################
# removed genomes
################################################################################

def test_removed_genomes_are_not_this_module_s_job(tmp_path):
    """
    gen-scg-hmms used to write its own `missed-accessions.tsv`. It now reports losses
    through the shared `write_removed_genomes_report`, so this program's account of
    what it lost has the same filename, columns, and wording as the main driver's and
    the search subcommands'. A second writer here would be a second answer.
    """
    assert not hasattr(outputs, "write_missed_accessions")
    assert not hasattr(outputs, "MISSED_FILENAME")


def test_removed_genomes_report_covers_a_gen_scg_hmms_run(tmp_path):
    from gtotree.utils.misc.summary_info import write_removed_genomes_report
    from gtotree.utils.misc.messaging import REMOVED_GENOMES_FILENAME
    from gtotree.utils.misc.stages import GenomeRemovalStage

    run_data = _run_data(tmp_path)
    run_data.ncbi_accs[0].mark_removed("not found in NCBI assembly data",
                                       GenomeRemovalStage.NCBI_LOOKUP)
    run_data.fasta_files[0].mark_removed("prodigal failed",
                                         GenomeRemovalStage.FASTA_PREP)

    write_removed_genomes_report(run_data)

    rows = _read(str(tmp_path / REMOVED_GENOMES_FILENAME))
    assert rows[0] == ["genome_id", "input", "source", "stage_removed",
                       "reason_removed"]
    assert rows[1] == ["G1", "-w Nitrospirota", "ncbi-accession (via -w)",
                       "ncbi-lookup", "not found in NCBI assembly data"]
    assert rows[2] == ["g2", "g2.fna", "nucleotide-fasta", "fasta-prep",
                       "prodigal failed"]


################################################################################
# target genomes
################################################################################

def _run_data(tmp_path):
    """A RunData shaped like a small gen-scg-hmms run: one `-w` accession, one fasta."""
    from gtotree.utils.misc.general import GenomeData, RunData

    run_data = RunData()
    run_data.run_files_dir = str(tmp_path)
    run_data.run_files_dir_rel = str(tmp_path)

    acc = GenomeData.from_acc("G1")
    acc.from_wanted_ref_tax = True
    acc.wanted_ref_tax_taxon = "Nitrospirota"
    acc.organism_name = "Nitrospira sp."
    run_data.ncbi_accs = [acc]

    run_data.fasta_files = [GenomeData.from_path("g2.fna", "nucleotide-fasta")]
    run_data.update_all_input_genomes()
    return run_data


def test_write_target_genomes(tmp_path):
    run_data = _run_data(tmp_path)
    path = outputs.write_target_genomes(str(tmp_path), ["G1", "g2"], run_data)

    rows = _read(path)
    # the same `input` and `source` columns the removed-genomes report uses, from the
    # same helpers, so the two tables can be read side by side
    assert rows[0] == ["genome_id", "input", "source", "organism_name"]
    assert rows[1] == ["G1", "-w Nitrospirota", "ncbi-accession (via -w)",
                       "Nitrospira sp."]
    assert rows[2] == ["g2", "g2.fna", "nucleotide-fasta", "NA"]


def test_write_target_genomes_defaults_to_na_for_an_unknown_genome(tmp_path):
    path = outputs.write_target_genomes(str(tmp_path), ["ghost"],
                                        _run_data(tmp_path))
    assert _read(path)[1] == ["ghost", "NA", "NA", "NA"]


################################################################################
# pfam version + naming
################################################################################

def test_write_pfam_version(tmp_path):
    path = outputs.write_pfam_version(str(tmp_path), "38.2")
    assert open(path).read().strip() == "38.2"


def test_default_hmm_filename_uses_named_output_dir(tmp_path):
    assert outputs.default_hmm_filename("My-SCG-Set", 42) == "My-SCG-Set.hmm"


def test_default_hmm_filename_falls_back_on_default_dir():
    assert outputs.default_hmm_filename("gtt-gen-scg-hmms-output", 42) == \
        "wanted-42-scg-targets.hmm"


def test_default_hmm_filename_strips_trailing_slash():
    assert outputs.default_hmm_filename("My-Set/", 3) == "My-Set.hmm"


################################################################################
# atomicity
################################################################################

def test_atomic_write_leaves_nothing_on_failure(tmp_path):
    path = tmp_path / "out.tsv"

    with pytest.raises(RuntimeError):
        with outputs._atomic_write(str(path)) as handle:
            handle.write("partial")
            raise RuntimeError("boom")

    assert not path.exists()
    assert not (tmp_path / "out.tsv.part").exists()


def test_atomic_write_moves_into_place_on_success(tmp_path):
    path = tmp_path / "out.tsv"
    with outputs._atomic_write(str(path)) as handle:
        handle.write("done")

    assert path.read_text() == "done"
    assert not (tmp_path / "out.tsv.part").exists()
