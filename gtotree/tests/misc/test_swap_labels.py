import os
import pytest  # type: ignore
from gtotree.utils.misc import seqs
from gtotree.utils.misc.general import RunData
from gtotree.utils.misc.seqs import swap_labels_in_alignment


ALIGNMENT = ">G0\nMMMMM\n>G1\nKKKKK\n"


def _run_data(tmp_path):
    out_dir = tmp_path / "output"
    run_files_dir = out_dir / "run-files"
    run_files_dir.mkdir(parents=True)

    rd = RunData()
    rd.output_dir = str(out_dir)
    rd.run_files_dir = str(run_files_dir)
    rd.general_ext = ".faa"
    rd.updating_headers = True
    rd.mapping_dict = {"G0": "Genus species G0", "G1": "Genus species G1"}
    rd.concatenated_alignment_path = str(out_dir / "aligned-SCGs.faa")
    return rd


def _paths(rd):
    return (os.path.join(rd.output_dir, "aligned-SCGs.faa"),
            os.path.join(rd.run_files_dir, "aligned-SCGs.faa"),
            os.path.join(rd.output_dir, "aligned-SCGs-mod-names.faa"))


class TestSwapLabels:

    def test_relabels_and_moves_the_original(self, tmp_path):
        rd = _run_data(tmp_path)
        orig, backup, new = _paths(rd)
        open(orig, "w").write(ALIGNMENT)

        rd = swap_labels_in_alignment(rd)

        assert rd.final_alignment_path == new
        assert ">Genus species G0" in open(new).read()
        # the original moves out of the output dir, it isn't left in both places
        assert os.path.exists(backup)
        assert not os.path.exists(orig)

    def test_falls_back_to_the_moved_original(self, tmp_path):
        """
        A run killed after the move but before headers_updated was persisted
        re-enters here with the original already in run-files/.
        """
        rd = _run_data(tmp_path)
        orig, backup, new = _paths(rd)
        open(backup, "w").write(ALIGNMENT)

        rd = swap_labels_in_alignment(rd)

        assert ">Genus species G1" in open(new).read()
        assert os.path.exists(backup)
        assert not os.path.exists(orig)

    def test_is_idempotent(self, tmp_path):
        rd = _run_data(tmp_path)
        orig, backup, new = _paths(rd)
        open(orig, "w").write(ALIGNMENT)

        rd = swap_labels_in_alignment(rd)
        first = open(new).read()
        rd = swap_labels_in_alignment(rd)

        assert open(new).read() == first
        assert os.path.exists(backup)
        assert not os.path.exists(orig)

    def test_an_empty_original_falls_through_to_the_backup(self, tmp_path):
        rd = _run_data(tmp_path)
        orig, backup, new = _paths(rd)
        open(orig, "w").close()
        open(backup, "w").write(ALIGNMENT)

        rd = swap_labels_in_alignment(rd)

        assert ">Genus species G0" in open(new).read()
        assert not os.path.exists(orig)

    def test_missing_from_both_places_raises(self, tmp_path):
        rd = _run_data(tmp_path)
        with pytest.raises(FileNotFoundError):
            swap_labels_in_alignment(rd)

    def test_an_interrupted_write_leaves_the_original_in_place(self, tmp_path,
                                                               monkeypatch):
        """
        The relabeled file is written before the move, so an interrupted write must
        leave the source exactly where a retry expects to find it.
        """
        rd = _run_data(tmp_path)
        orig, backup, new = _paths(rd)
        open(orig, "w").write(ALIGNMENT)

        def exploding_parse(*a, **kw):
            raise RuntimeError("killed mid-write")

        monkeypatch.setattr(seqs.SeqIO, "parse", exploding_parse)
        with pytest.raises(RuntimeError):
            swap_labels_in_alignment(rd)
        monkeypatch.undo()

        assert not os.path.exists(new)
        assert not os.path.exists(new + ".part")
        assert os.path.exists(orig)
        assert not os.path.exists(backup)

        # and the retry succeeds
        rd = swap_labels_in_alignment(rd)
        assert ">Genus species G0" in open(new).read()
