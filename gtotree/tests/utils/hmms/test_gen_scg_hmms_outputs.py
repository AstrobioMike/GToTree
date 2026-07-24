import os
import pytest # type: ignore
from gtotree.utils.hmms import gen_scg_hmms_outputs as outputs
from gtotree.utils.hmms.gen_scg_hmms import PfamProfileInfo


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
# missed accessions
################################################################################

def test_write_missed_accessions(tmp_path):
    path = outputs.write_missed_accessions(
        str(tmp_path), [("GCA_1", "not found"), ("GCA_2", "download failed")])

    rows = _read(path)
    assert rows[0] == ["accession", "why_missing"]
    assert rows[1] == ["GCA_1", "not found"]
    assert rows[2] == ["GCA_2", "download failed"]


def test_no_missed_file_when_nothing_missed(tmp_path):
    """No file at all is clearer than an empty one with just a header."""
    assert outputs.write_missed_accessions(str(tmp_path), []) is None
    assert not (tmp_path / outputs.MISSED_FILENAME).exists()


################################################################################
# target genomes
################################################################################

def test_write_target_genomes(tmp_path):
    path = outputs.write_target_genomes(
        str(tmp_path), ["G1"], {"G1": "GTDB:Nitrospirota"}, {"G1": "Nitrospira sp."})

    rows = _read(path)
    assert rows[0] == ["accession", "source", "organism_name"]
    assert rows[1] == ["G1", "GTDB:Nitrospirota", "Nitrospira sp."]


def test_write_target_genomes_defaults_to_na(tmp_path):
    path = outputs.write_target_genomes(str(tmp_path), ["G1"])
    assert _read(path)[1] == ["G1", "NA", "NA"]


def test_write_target_genomes_handles_none_organism(tmp_path):
    path = outputs.write_target_genomes(str(tmp_path), ["G1"], {"G1": "fasta"},
                                        {"G1": None})
    assert _read(path)[1] == ["G1", "fasta", "NA"]


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
