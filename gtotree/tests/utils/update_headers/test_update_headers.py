"""
Tests for `gtt update-headers`.

The behavior worth pinning down is that this re-labels a finished run *correctly*
without touching it, in the three situations that differ materially: a run that never
re-labelled (outputs keyed by genome id), one that did (outputs keyed by the old
labels, and the original-ID alignment moved to run-files/), and one whose swap failed
(mapping dict populated, but outputs still keyed by genome id).

Everything taxonomy-related is exercised through a stubbed handler rather than the real
Parquet assets, which aren't available in the test environment -- the handlers
themselves are covered by their own tests, and what matters here is that this program
hands them a *cleared* mapping dict.
"""

import os
import types
import pytest  # type: ignore

from gtotree.tests.utils.update_headers.conftest import (
    ACCESSIONS, DROPPED_ACCESSION, build_completed_run)
from gtotree.utils.update_headers import update_headers_run as run_lib
from gtotree.utils.update_headers import update_headers_cli as cli
from gtotree.utils.update_headers.update_headers_run import UpdateHeadersError


def _args(input_dir, output_dir, **overrides):
    args = types.SimpleNamespace(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        output_prefix="",
        mapping_file=None,
        add_gtdb_tax=False,
        add_ncbi_tax=False,
        lineage=cli.DEFAULT_LINEAGE,
        force_overwrite=False,
    )
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def _fasta_headers(path):
    return [line[1:].strip() for line in open(path) if line.startswith(">")]


def _tsv(path):
    lines = open(path).read().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def _outputs(tmp_path, input_dir, **overrides):
    """Run the whole thing and return (out_dir, args, run_data)."""
    out_dir = overrides.pop("output_dir", tmp_path / "relabeled")
    args = _args(input_dir, out_dir, **overrides)
    args = cli.check_args(args)

    run_data = run_lib.load_completed_run(args.input_dir)
    run_data = run_lib.rehome_run_data(run_data, args.input_dir)

    previous = run_lib.labels_in_completed_outputs(run_data)

    args, run_data = cli.build_new_mapping_dict(args, run_data)
    new_labels = run_lib.new_labels_for(run_data)

    written_dir, _written, _total, _changed = cli.write_outputs(
        args, run_data, previous, new_labels)

    return written_dir, args, run_data


################################################################################
# reading the completed run
################################################################################

class TestLoadingACompletedRun:

    def test_loads_a_finished_run(self, completed_run):
        run_data = run_lib.load_completed_run(completed_run)
        assert [gd.id for gd in run_data.ncbi_accs] == ACCESSIONS + [DROPPED_ACCESSION]

    def test_a_missing_directory_says_so(self, tmp_path):
        with pytest.raises(UpdateHeadersError, match="doesn't exist"):
            run_lib.load_completed_run(str(tmp_path / "nope"))

    def test_a_directory_that_isnt_a_gtotree_run_says_so(self, tmp_path):
        plain = tmp_path / "just-a-dir"
        plain.mkdir()
        with pytest.raises(UpdateHeadersError, match="doesn't look like a GToTree"):
            run_lib.load_completed_run(str(plain))

    def test_a_run_without_an_alignment_stage_is_refused(self, tmp_path):
        """
        An interrupted run has a run-data.json but nothing to re-label, and pointing
        the user back at `-R` is a more useful answer than a missing-file error later.
        """
        from gtotree.utils.misc.general import read_run_data, write_run_data
        from gtotree.utils.misc.stages import PipelineStage

        out_dir = build_completed_run(tmp_path)
        run_data = read_run_data(os.path.join(out_dir, "run-files", "run-data.json"))
        run_data.completed_stages.pop(PipelineStage.CONCATENATE_SCG_SETS)
        write_run_data(run_data)

        with pytest.raises(UpdateHeadersError, match="concatenated alignment"):
            run_lib.load_completed_run(out_dir)

    def test_corrupt_run_data_is_a_clear_error(self, tmp_path):
        out_dir = build_completed_run(tmp_path)
        with open(os.path.join(out_dir, "run-files", "run-data.json"), "w") as f:
            f.write("{not json")

        with pytest.raises(UpdateHeadersError, match="interrupted while saving"):
            run_lib.load_completed_run(out_dir)


class TestRehoming:
    """
    Every path in run-data.json is absolute, so a run that's been moved or renamed is
    the normal case for a tool people reach for weeks after the fact.
    """

    def test_a_moved_run_still_resolves_its_outputs(self, tmp_path):
        original = build_completed_run(tmp_path / "first")
        moved = tmp_path / "second" / "gtotree-output"
        os.makedirs(tmp_path / "second")
        os.rename(original, moved)

        run_data = run_lib.load_completed_run(str(moved))
        run_data = run_lib.rehome_run_data(run_data, str(moved))

        assert run_data.run_files_dir == os.path.join(str(moved), "run-files")
        assert os.path.isfile(run_data.final_tree_path)
        assert os.path.isfile(run_lib.original_id_alignment_path(run_data))

    def test_a_renamed_run_still_finds_its_tree(self, tmp_path):
        """
        The tree is named after the output directory, so a rename breaks the basename
        lookup and only the single-.tre fallback can find it.
        """
        original = build_completed_run(tmp_path, name="gtotree-output")
        renamed = str(tmp_path / "renamed-output")
        os.rename(original, renamed)

        run_data = run_lib.load_completed_run(renamed)
        run_data = run_lib.rehome_run_data(run_data, renamed)

        assert os.path.basename(run_data.final_tree_path) == "gtotree-output.tre"

    def test_the_dangling_working_dir_paths_are_dropped(self, completed_run):
        """
        A finished run deletes its working directory, so the NCBI sub-table recorded
        there is gone; carrying the stale path forward would blow up the summary table.
        """
        run_data = run_lib.load_completed_run(completed_run)
        assert run_data.ncbi_sub_table_path  # recorded by the original run

        run_data = run_lib.rehome_run_data(run_data, completed_run)
        assert run_data.ncbi_sub_table_path == ""


################################################################################
# labels
################################################################################

class TestPreviousLabels:

    def test_a_run_that_never_relabeled_is_keyed_by_genome_id(self, completed_run):
        run_data = run_lib.load_completed_run(completed_run)
        previous = run_lib.labels_in_completed_outputs(run_data)
        assert previous[ACCESSIONS[0]] == ACCESSIONS[0]

    def test_a_relabeled_run_reports_its_labels(self, relabeled_run):
        run_data = run_lib.load_completed_run(relabeled_run)
        previous = run_lib.labels_in_completed_outputs(run_data)
        assert previous["GCF_000005845.2"] == "E_coli_K12"

    def test_a_failed_swap_still_reports_genome_ids(self, tmp_path):
        """
        The mapping dict is populated but the outputs kept the original IDs, so the
        stage's completion is what settles this, not the presence of a mapping.
        """
        out_dir = build_completed_run(tmp_path, swapped=False, mapping={
            "GCF_000005845.2": "E_coli_K12"})

        run_data = run_lib.load_completed_run(out_dir)
        previous = run_lib.labels_in_completed_outputs(run_data)
        assert previous["GCF_000005845.2"] == "GCF_000005845.2"


class TestBuildingTheNewMapping:

    def test_the_previous_mapping_is_cleared_first(self, relabeled_run, monkeypatch):
        """
        Both taxonomy handlers skip genomes already in the mapping dict. Carrying the
        finished run's mapping in would mean a `-D` re-label silently did nothing,
        which is the whole failure this program exists to avoid.
        """
        seen = {}

        def fake_gtdb(args, run_data):
            seen["dict_at_entry"] = dict(run_data.mapping_dict)
            for acc in ACCESSIONS:
                run_data.mapping_dict[acc] = f"{acc}_Bacteria_Pseudomonadota"
            return run_data

        monkeypatch.setattr(
            "gtotree.utils.gtdb.handle_gtdb_tax_info.update_mapping_dict_with_gtdb_tax_info",
            fake_gtdb)

        run_data = run_lib.load_completed_run(relabeled_run)
        args = _args(relabeled_run, "unused", add_gtdb_tax=True)
        _args_out, run_data = cli.build_new_mapping_dict(args, run_data)

        assert seen["dict_at_entry"] == {}
        assert run_data.mapping_dict["GCF_000005845.2"].endswith("_Pseudomonadota")

    def test_no_label_flags_falls_back_to_genome_ids(self, relabeled_run):
        """
        Running with nothing specified is a legitimate way to put the original IDs
        back on a run that was labeled.
        """
        run_data = run_lib.load_completed_run(relabeled_run)
        args = _args(relabeled_run, "unused")
        _args_out, run_data = cli.build_new_mapping_dict(args, run_data)

        new_labels = run_lib.new_labels_for(run_data)
        assert new_labels["GCF_000005845.2"] == "GCF_000005845.2"


################################################################################
# outputs
################################################################################

class TestWrittenOutputs:

    def test_the_alignment_carries_the_new_labels(self, tmp_path, relabeled_run):
        mapping = tmp_path / "new-names.tsv"
        mapping.write_text("GCF_000005845.2\tK12\nGCA_000009065.1\tSakai\n")

        out_dir, _args_out, _rd = _outputs(tmp_path, relabeled_run,
                                           mapping_file=str(mapping))

        headers = _fasta_headers(os.path.join(out_dir,
                                              "aligned-SCGs-mod-names.faa"))
        assert headers == ["K12", "Sakai", "GCF_000006945.2"]

    def test_the_tree_carries_the_new_labels(self, tmp_path, relabeled_run):
        """
        The tree can only be re-labeled through the previous labels -- there is no
        original-ID tree, since the tree is built from whatever the alignment carried.
        """
        mapping = tmp_path / "new-names.tsv"
        mapping.write_text("GCF_000005845.2\tK12\nGCA_000009065.1\tSakai\n")

        out_dir, _args_out, _rd = _outputs(tmp_path, relabeled_run,
                                           mapping_file=str(mapping))

        newick = open(os.path.join(out_dir, "gtotree-output.tre")).read()
        assert "K12" in newick and "Sakai" in newick
        assert "E_coli_K12" not in newick
        # the untouched genome keeps its label, and support values survive
        assert "GCF_000006945.2" in newick
        assert "0.98" in newick

    def test_an_unlabeled_run_relabels_from_genome_ids(self, tmp_path, completed_run):
        mapping = tmp_path / "new-names.tsv"
        mapping.write_text("GCF_000005845.2\tK12\n")

        out_dir, _args_out, _rd = _outputs(tmp_path, completed_run,
                                           mapping_file=str(mapping))

        newick = open(os.path.join(out_dir, "gtotree-output.tre")).read()
        assert "K12" in newick
        assert "GCF_000005845.2" not in newick

    def test_the_summary_table_holds_the_new_labels_and_the_taxids(
            self, tmp_path, relabeled_run):
        mapping = tmp_path / "new-names.tsv"
        mapping.write_text("GCF_000005845.2\tK12\n")

        out_dir, _args_out, _rd = _outputs(tmp_path, relabeled_run,
                                           mapping_file=str(mapping))

        rows = {r["genome_id"]: r
                for r in _tsv(os.path.join(out_dir, "genomes-summary-info.tsv"))}
        assert rows["GCF_000005845.2"]["label"] == "K12"
        # the taxid rode along on the GenomeData; the sub-table it originally came
        # from was deleted with the working directory when the run finished
        assert rows["GCF_000005845.2"]["taxid"] == "511145"
        assert rows[DROPPED_ACCESSION]["in_final_tree"] == "No"

    def test_the_label_map_records_both_sides(self, tmp_path, relabeled_run):
        mapping = tmp_path / "new-names.tsv"
        mapping.write_text("GCF_000005845.2\tK12\n")

        out_dir, _args_out, _rd = _outputs(tmp_path, relabeled_run,
                                           mapping_file=str(mapping))

        rows = {r["genome_id"]: r
                for r in _tsv(os.path.join(out_dir, "label-map.tsv"))}
        assert rows["GCF_000005845.2"]["previous_label"] == "E_coli_K12"
        assert rows["GCF_000005845.2"]["new_label"] == "K12"
        # a genome the original run dropped is still accounted for
        assert rows[DROPPED_ACCESSION]["new_label"] == DROPPED_ACCESSION

    def test_the_prefix_is_applied_to_every_output(self, tmp_path, completed_run):
        out_dir, _args_out, _rd = _outputs(tmp_path, completed_run,
                                           output_prefix="family-level")

        names = sorted(os.listdir(out_dir))
        assert names == ["family-level-aligned-SCGs-mod-names.faa",
                         "family-level-genomes-summary-info.tsv",
                         "family-level-gtotree-output.tre",
                         "family-level-label-map.tsv"]

    def test_nucleotide_mode_uses_the_right_extension(self, tmp_path):
        source = build_completed_run(tmp_path, nucleotide_mode=True)
        out_dir, _args_out, _rd = _outputs(tmp_path, source)

        assert os.path.isfile(os.path.join(out_dir, "aligned-SCGs-mod-names.fasta"))

    def test_a_run_with_no_tree_still_produces_the_rest(self, tmp_path):
        source = build_completed_run(tmp_path, with_tree=False)
        out_dir, _args_out, _rd = _outputs(tmp_path, source)

        assert not [n for n in os.listdir(out_dir) if n.endswith(".tre")]
        assert os.path.isfile(os.path.join(out_dir, "aligned-SCGs-mod-names.faa"))

    def test_the_original_run_is_left_alone(self, tmp_path, relabeled_run):
        before = _snapshot(relabeled_run)

        mapping = tmp_path / "new-names.tsv"
        mapping.write_text("GCF_000005845.2\tK12\n")
        _outputs(tmp_path, relabeled_run, mapping_file=str(mapping))

        assert _snapshot(relabeled_run) == before

    def test_missing_run_files_falls_back_to_the_final_alignment(self, tmp_path):
        """
        People delete run-files/ to save space. The finished alignment can still be
        re-labeled through the previous labels, the same way the tree is.
        """
        import shutil

        source = build_completed_run(tmp_path, mapping={
            "GCF_000005845.2": "E_coli_K12",
            "GCA_000009065.1": "E_coli_Sakai",
            "GCF_000006945.2": "S_enterica_LT2"})

        run_data = run_lib.load_completed_run(source)
        run_data = run_lib.rehome_run_data(run_data, source)
        previous = run_lib.labels_in_completed_outputs(run_data)

        shutil.rmtree(os.path.join(source, "run-files"))

        out_path = tmp_path / "out.faa"
        run_lib.write_updated_alignment(
            run_data, previous,
            {acc: f"new-{i}" for i, acc in enumerate(ACCESSIONS)},
            str(out_path))

        assert _fasta_headers(out_path) == ["new-0", "new-1", "new-2"]


def _snapshot(directory):
    """Every file under `directory`, with its contents, for an untouched-check."""
    snapshot = {}
    for root, _dirs, files in os.walk(directory):
        for name in files:
            path = os.path.join(root, name)
            snapshot[os.path.relpath(path, directory)] = open(path, "rb").read()
    return snapshot


################################################################################
# argument validation
################################################################################

class TestArgChecks:

    def test_both_taxonomies_at_once_is_refused(self, tmp_path):
        args = _args(tmp_path, tmp_path / "out",
                     add_gtdb_tax=True, add_ncbi_tax=True)
        with pytest.raises(UpdateHeadersError, match="one or the other"):
            cli.check_args(args)

    def test_an_unknown_rank_is_refused(self, tmp_path):
        args = _args(tmp_path, tmp_path / "out",
                     add_gtdb_tax=True, lineage="domain,subphylum")
        with pytest.raises(UpdateHeadersError, match="not an accepted taxonomic rank"):
            cli.check_args(args)

    def test_a_custom_lineage_without_a_taxonomy_flag_is_refused(self, tmp_path):
        args = _args(tmp_path, tmp_path / "out", lineage="domain,family")
        with pytest.raises(UpdateHeadersError, match="which taxonomy"):
            cli.check_args(args)

    def test_writing_into_the_run_being_read_is_refused(self, tmp_path):
        """
        The one thing this program promises is that it doesn't touch the run it reads.
        """
        args = _args(tmp_path, str(tmp_path) + "/")
        with pytest.raises(UpdateHeadersError, match="never writes into the run"):
            cli.check_args(args)

    def test_an_existing_output_dir_is_refused_without_force(self, tmp_path,
                                                             completed_run):
        existing = tmp_path / "already-here"
        existing.mkdir()

        args = _args(completed_run, existing)
        from gtotree.utils.misc.general import OutputDirExistsError
        with pytest.raises(OutputDirExistsError):
            cli.setup_output_dir(args)

    def test_force_overwrite_replaces_it(self, tmp_path, completed_run):
        existing = tmp_path / "already-here"
        existing.mkdir()
        (existing / "stale.txt").write_text("old")

        args = _args(completed_run, existing, force_overwrite=True)
        cli.setup_output_dir(args)

        assert not (existing / "stale.txt").exists()


################################################################################
# iToL regeneration
################################################################################

class TestItolRegeneration:

    def test_itol_files_are_rebuilt_with_the_new_labels(self, tmp_path):
        source = build_completed_run(tmp_path, target_pfams_file="targets.txt")

        pfam_dir = os.path.join(source, "pfam-search-results")
        os.makedirs(pfam_dir)
        with open(os.path.join(pfam_dir, "pfam-hit-counts.tsv"), "w") as f:
            f.write("genome_id\ttotal_gene_count\tPF00789\n")
            f.write(f"{ACCESSIONS[0]}\t4000\t2\n")
            f.write(f"{ACCESSIONS[1]}\t4100\t0\n")

        mapping = tmp_path / "new-names.tsv"
        mapping.write_text(f"{ACCESSIONS[0]}\tK12\n")

        out_dir, _args_out, _rd = _outputs(
            tmp_path, source, mapping_file=str(mapping))

        itol_file = os.path.join(out_dir, "pfam-search-results", "iToL-files",
                                 "PF00789-iToL.txt")
        assert "K12" in open(itol_file).read()

    def test_nothing_is_written_when_the_run_had_no_searches(self, tmp_path,
                                                             completed_run):
        out_dir, _args_out, _rd = _outputs(tmp_path, completed_run)
        assert not os.path.isdir(os.path.join(out_dir, "pfam-search-results"))
