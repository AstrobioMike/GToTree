"""
Tests for the setup layer of `gtt search-annotations`.

The two things worth pinning down here are the ones that would be easy to get subtly
wrong and hard to notice:

  * the dependency checks are *narrower* than the main driver's, deliberately -- a
    regression that reintroduced the muscle/trimal requirement would break working
    setups for no reason
  * the flattened output layout, since the search helpers all write to paths derived
    from run_data attributes, and getting one wrong scatters files
"""

import os

import pytest  # type: ignore

from gtotree.utils.target_search.target_search_setup import (
    TargetSearchError,
    check_args,
    check_dependencies,
    check_env_vars,
    setup_output_dir,
    build_run_data,
    ensure_processing_dirs,
    ensure_reference_data,
    validate_input_files,
    make_args,
    WORKING_DIR_NAME,
    RUN_DATA_FILENAME,
)
from gtotree.utils.target_search.target_search_spec import get_spec
from gtotree.utils.misc.general import OutputDirExistsError, adopt_genome_progress


@pytest.fixture
def spec():
    return get_spec("pfam")


################################################################################
# check_args
################################################################################

def test_requires_at_least_one_genome_source(spec):
    args = make_args(target_pfams_file="t.txt")
    with pytest.raises(TargetSearchError, match="input genomes"):
        check_args(args, spec)


def test_requires_a_targets_file(spec):
    args = make_args(amino_acid_files="aa.txt")
    with pytest.raises(TargetSearchError, match="which Pfams"):
        check_args(args, spec)


def test_ko_spec_names_kos_in_the_missing_targets_message():
    """The message has to name the right target type, since the spec drives it."""
    args = make_args(amino_acid_files="aa.txt")
    with pytest.raises(TargetSearchError, match="which KOs"):
        check_args(args, get_spec("ko"))


@pytest.mark.parametrize("source_flag", ["ncbi_accessions", "genbank_files",
                                         "fasta_files", "amino_acid_files",
                                         "wanted_ref_tax"])
def test_any_single_genome_source_is_enough(spec, source_flag):
    args = make_args(target_pfams_file="t.txt", **{source_flag: "x"})
    assert check_args(args, spec) is args


def test_rejects_resume_with_force_overwrite(spec):
    args = make_args(amino_acid_files="aa.txt", target_pfams_file="t.txt",
                     resume=True, force_overwrite=True)
    with pytest.raises(TargetSearchError, match="can't be used together"):
        check_args(args, spec)


def test_rejects_zero_jobs(spec):
    args = make_args(amino_acid_files="aa.txt", target_pfams_file="t.txt", num_jobs=0)
    with pytest.raises(TargetSearchError, match="--num-jobs"):
        check_args(args, spec)


def test_rejects_dangling_target_rank(spec):
    args = make_args(amino_acid_files="aa.txt", target_pfams_file="t.txt",
                     target_rank="genus")
    with pytest.raises(TargetSearchError, match="--target-rank"):
        check_args(args, spec)


def test_rejects_dangling_derep_rank(spec):
    args = make_args(amino_acid_files="aa.txt", target_pfams_file="t.txt",
                     derep_rank="genus")
    with pytest.raises(TargetSearchError, match="--derep-rank"):
        check_args(args, spec)


def test_default_derep_rank_is_not_treated_as_dangling(spec):
    """`--derep-rank` has a non-None default, so only a user-changed value is dangling."""
    args = make_args(amino_acid_files="aa.txt", target_pfams_file="t.txt",
                     derep_rank="auto")
    assert check_args(args, spec) is args


def test_no_minimum_genome_count(spec, tmp_path, in_tmp_cwd, listing, write_genome):
    """
    Searching a single genome for one Pfam is legitimate, so the main driver's
    four-genome floor must not apply here.
    """
    genome = write_genome("only", ["PF90001"])
    args = make_args(amino_acid_files=listing("aa.txt", [genome]),
                     target_pfams_file=listing("t.txt", ["PF90001"]))
    check_args(args, spec)
    args = validate_input_files(args, spec)
    out_dir, work_dir = setup_output_dir(args, spec)
    run_data = build_run_data(args, spec, out_dir, work_dir)
    assert len(run_data.all_input_genomes) == 1


################################################################################
# dependency checks
################################################################################

def test_does_not_require_muscle_or_trimal(spec, monkeypatch):
    """
    The whole reason these subcommands don't reuse `preflight_checks` -- neither aligns
    nor trims anything, so demanding the aligner would fail working setups.
    """
    seen = []

    def fake_which(cmd):
        seen.append(cmd)
        return "/usr/bin/" + cmd

    monkeypatch.setattr("shutil.which", fake_which)
    check_dependencies(make_args(amino_acid_files="aa.txt"), spec)

    assert "muscle" not in seen
    assert "trimal" not in seen


def test_amino_acid_only_run_does_not_need_prodigal(spec, monkeypatch):
    monkeypatch.setattr("shutil.which", lambda cmd: None)
    # would raise if prodigal were required
    check_dependencies(make_args(amino_acid_files="aa.txt"), spec)


@pytest.mark.parametrize("source_flag", ["fasta_files", "genbank_files",
                                         "ncbi_accessions", "wanted_ref_tax"])
def test_prodigal_required_when_gene_calling_is_possible(spec, monkeypatch, source_flag):
    """
    Nucleotide fastas always call genes; NCBI assemblies with no protein file and
    GenBank files with no usable CDS translations both fall back to prodigal.
    """
    monkeypatch.setattr("shutil.which", lambda cmd: None)
    args = make_args(**{source_flag: "x"})
    with pytest.raises(TargetSearchError, match="prodigal"):
        check_dependencies(args, spec)


def test_ko_run_requires_kofamscan(monkeypatch):
    """The KO search shells out to exec_annotation, so it's a hard requirement."""
    monkeypatch.setattr("shutil.which", lambda cmd: None)
    with pytest.raises(TargetSearchError, match="exec_annotation"):
        check_dependencies(make_args(amino_acid_files="aa.txt"), get_spec("ko"))


def test_env_var_check_catches_an_unset_data_dir(spec, monkeypatch):
    """
    The search helpers read these straight out of os.environ; without this check the
    failure surfaces as a bare KeyError inside a worker thread.
    """
    monkeypatch.delenv("Pfam_data_dir", raising=False)
    with pytest.raises(TargetSearchError, match="Pfam_data_dir"):
        check_env_vars(spec)


def test_env_var_check_passes_when_set(spec, monkeypatch):
    monkeypatch.setenv("Pfam_data_dir", "/somewhere")
    check_env_vars(spec)


################################################################################
# output directory
################################################################################

def test_creates_the_flattened_result_dirs(spec, tmp_path):
    args = make_args(output_dir=str(tmp_path / "out"))
    out_dir, work_dir = setup_output_dir(args, spec)

    assert work_dir == os.path.join(out_dir, WORKING_DIR_NAME)
    for sub in spec.result_subdirs:
        assert os.path.isdir(os.path.join(out_dir, sub))
    # the level a full GToTree run would have added
    assert not os.path.exists(os.path.join(out_dir, "pfam-search-results"))


def test_refuses_an_existing_output_dir(spec, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    args = make_args(output_dir=str(out))
    with pytest.raises(OutputDirExistsError, match="already exists"):
        setup_output_dir(args, spec)


def test_force_overwrite_clears_the_previous_run(spec, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "stale.txt").write_text("old")
    args = make_args(output_dir=str(out), force_overwrite=True)
    setup_output_dir(args, spec)
    assert not (out / "stale.txt").exists()


def test_resume_keeps_the_previous_run(spec, tmp_path):
    out = tmp_path / "out"
    (out / WORKING_DIR_NAME).mkdir(parents=True)
    (out / "keep.txt").write_text("mine")
    args = make_args(output_dir=str(out), resume=True)
    setup_output_dir(args, spec)
    assert (out / "keep.txt").exists()


def test_resume_into_a_fresh_directory_just_starts_over(spec, tmp_path):
    """A friendly no-op rather than a refusal -- there's nothing to mix up."""
    args = make_args(output_dir=str(tmp_path / "nope"), resume=True)
    out_dir, work_dir = setup_output_dir(args, spec)
    assert os.path.isdir(work_dir)


################################################################################
# run_data layout
################################################################################

@pytest.fixture
def built(spec, tmp_path, listing, write_genome):
    genomes = [write_genome(f"g{i}", ["PF90001"]) for i in (1, 2)]
    args = make_args(amino_acid_files=listing("aa.txt", genomes),
                     target_pfams_file=listing("t.txt", ["PF90001", "PF90002"]),
                     output_dir=str(tmp_path / "out"))
    args = validate_input_files(args, spec)
    out_dir, work_dir = setup_output_dir(args, spec)
    return args, build_run_data(args, spec, out_dir, work_dir), out_dir, work_dir


def test_results_dir_is_the_output_dir_itself(built, spec):
    _args, run_data, out_dir, _work = built
    assert spec.results_dir(run_data) == os.path.abspath(out_dir)


def test_run_files_dir_is_also_the_output_dir(built):
    """
    So run-level files like `removed-genomes.tsv` (written to run_files_dir by the
    shared helper) land at the top level rather than in a run-files subdirectory
    nothing else uses. (Per-type files like failed-<type>-targets.txt go to the type's
    results dir instead.)
    """
    _args, run_data, out_dir, _work = built
    assert run_data.run_files_dir == os.path.abspath(out_dir)


def test_intermediates_live_in_the_working_dir(built):
    _args, run_data, _out, work_dir = built
    assert run_data.tmp_dir.startswith(work_dir)
    assert run_data.ready_genome_files_dir.startswith(work_dir)
    assert run_data.run_data_path == os.path.join(work_dir, RUN_DATA_FILENAME)


def test_input_genomes_are_parsed_with_the_shared_factories(built):
    _args, run_data, _out, _work = built
    assert [gd.id for gd in run_data.all_input_genomes] == ["g1", "g2"]
    assert all(gd.source == "amino-acid-fasta" for gd in run_data.amino_acid_files)


def test_requested_target_count_is_recorded(built, spec):
    """Drives the "N of M requested targets found" line."""
    _args, run_data, _out, _work = built
    assert spec.total_targets(run_data) == 2


def test_nucleotide_mode_is_off(built):
    """There's no tree here, so there's nothing nucleotide mode could feed."""
    _args, run_data, _out, _work = built
    assert run_data.nucleotide_mode is False
    assert run_data.general_ext == ".faa"


def test_only_the_used_processing_dirs_are_created(built, spec):
    _args, run_data, _out, _work = built
    assert os.path.isdir(run_data.AA_processing_dir)
    assert not run_data.genbank_processing_dir
    assert not run_data.fasta_processing_dir


def test_resumed_build_reparses_the_genome_set(built, spec, listing, write_genome):
    """
    Regression test. `build_run_data` used to adopt the previous RunData wholesale on
    resume, which meant edits to the input listing were invisible: a genome added
    between runs was never processed, AND the fingerprint -- built from run_data --
    compared the old set against itself and saw nothing to refuse over. The set must
    always come from the current input files.
    """
    args, run_data, out_dir, work_dir = built

    extra = write_genome("g3", ["PF90001"])
    paths = [gd.provided_path for gd in run_data.all_input_genomes] + [str(extra)]
    args.amino_acid_files = listing("aa.txt", paths)

    rebuilt = build_run_data(args, spec, out_dir, work_dir, previous=run_data)
    assert [gd.id for gd in rebuilt.all_input_genomes] == ["g1", "g2", "g3"]


def test_resumed_build_carries_run_level_target_state(built, spec):
    """
    What `previous` legitimately contributes: the resolved target set, so the
    multi-minute Pfam collection doesn't re-run.
    """
    args, run_data, out_dir, work_dir = built
    run_data.found_pfam_targets = ["PF90001"]
    run_data.additional_pfam_searching_done = True

    rebuilt = build_run_data(args, spec, out_dir, work_dir, previous=run_data)
    assert rebuilt.found_pfam_targets == ["PF90001"]
    assert rebuilt.additional_pfam_searching_done is True


def test_adopt_genome_progress_carries_flags_by_id(built, spec):
    """
    The mechanism behind genome-level resume: flags move across by genome ID, and
    `genomes_needing_processing` filters on exactly those.
    """
    args, run_data, out_dir, work_dir = built
    run_data.all_input_genomes[0].mark_processing_done()
    run_data.all_input_genomes[0].mark_pfam_search_done()
    run_data.all_input_genomes[0].num_genes = 42

    rebuilt = build_run_data(args, spec, out_dir, work_dir, previous=run_data)
    carried = adopt_genome_progress(rebuilt, run_data)

    assert carried == 2
    assert rebuilt.all_input_genomes[0].pfam_search_done is True
    assert rebuilt.all_input_genomes[0].num_genes == 42
    assert rebuilt.all_input_genomes[1].pfam_search_done is False


def test_adopt_genome_progress_leaves_new_genomes_untouched(built, spec, listing,
                                                            write_genome):
    """A genome that wasn't in the previous run keeps its fresh flags and gets done."""
    args, run_data, out_dir, work_dir = built
    for gd in run_data.all_input_genomes:
        gd.mark_processing_done()
        gd.mark_pfam_search_done()

    extra = write_genome("g3", ["PF90001"])
    paths = [gd.provided_path for gd in run_data.all_input_genomes] + [str(extra)]
    args.amino_acid_files = listing("aa.txt", paths)

    rebuilt = build_run_data(args, spec, out_dir, work_dir, previous=run_data)
    adopt_genome_progress(rebuilt, run_data)

    by_id = {gd.id: gd for gd in rebuilt.all_input_genomes}
    assert by_id["g1"].pfam_search_done is True
    assert by_id["g3"].pfam_search_done is False


def test_adopt_genome_progress_never_overwrites_identity(built, spec):
    """
    Progress moves across; where the genome came from does not. Otherwise a stale
    `full_path` from the previous run could outlive the file it names.
    """
    args, run_data, out_dir, work_dir = built
    run_data.all_input_genomes[0].full_path = "/gone/stale.faa"

    rebuilt = build_run_data(args, spec, out_dir, work_dir, previous=run_data)
    adopt_genome_progress(rebuilt, run_data)

    assert rebuilt.all_input_genomes[0].full_path != "/gone/stale.faa"


def test_adopt_genome_progress_skips_a_changed_source(built, spec):
    """
    Same basename arriving as a different input type is a different genome; reusing its
    flags would skip work that genuinely needs redoing.
    """
    args, run_data, out_dir, work_dir = built
    run_data.all_input_genomes[0].source = "genbank-file"
    run_data.all_input_genomes[0].mark_pfam_search_done()

    rebuilt = build_run_data(args, spec, out_dir, work_dir, previous=run_data)
    adopt_genome_progress(rebuilt, run_data)

    assert rebuilt.all_input_genomes[0].pfam_search_done is False


def test_adopt_genome_progress_with_no_previous_run_is_a_noop(built, spec):
    _args, run_data, _out, _work = built
    assert adopt_genome_progress(run_data, None) == 0


################################################################################
# processing directories
################################################################################

"""
Regression tests for a real download failure.

`--wanted-ref-tax` accessions are merged into `run_data.ncbi_accs` during phase 1,
which is AFTER `build_run_data` runs. Directory creation used to be gated on
`run_data.ncbi_accs` being non-empty at build time, so a run whose only NCBI accessions
came from `-w` left `ncbi_processing_dir` as "". `prepare_accession` concatenates that
into its download path, so every genome was downloaded to "/<acc>_protein.faa" at the
filesystem root, failed on permissions, and got swallowed by that function's bare
`except:` -- surfacing only as "acc download failed" for every genome.
"""


def test_wanted_ref_tax_only_run_gets_a_processing_dir(spec, tmp_path, listing,
                                                       write_genome):
    args = make_args(wanted_ref_tax="Nitrospirota",
                     target_pfams_file=listing("t.txt", ["PF90001"]),
                     output_dir=str(tmp_path / "out"))
    args = validate_input_files(args, spec)
    out_dir, work_dir = setup_output_dir(args, spec)
    run_data = build_run_data(args, spec, out_dir, work_dir)

    # nothing to download yet, so nothing is set up yet
    assert run_data.ncbi_processing_dir == ""

    # ... this is what phase 1 does
    run_data.merge_wanted_ref_tax_accessions(["GCF_000008865.2", "GCF_000005845.2"])
    ensure_processing_dirs(run_data)

    assert run_data.ncbi_processing_dir
    assert os.path.isdir(run_data.ncbi_processing_dir)
    # the actual failure mode: a download target at the filesystem root
    assert not run_data.ncbi_processing_dir.startswith("/GCF")


def test_ensure_processing_dirs_is_idempotent(spec, tmp_path, listing, write_genome):
    """Called twice on every run -- once at build, once after the genome set is final."""
    genome = write_genome("g1", ["PF90001"])
    args = make_args(amino_acid_files=listing("aa.txt", [genome]),
                     target_pfams_file=listing("t.txt", ["PF90001"]),
                     output_dir=str(tmp_path / "out"))
    args = validate_input_files(args, spec)
    out_dir, work_dir = setup_output_dir(args, spec)
    run_data = build_run_data(args, spec, out_dir, work_dir)

    first = run_data.AA_processing_dir
    ensure_processing_dirs(run_data)
    ensure_processing_dirs(run_data)

    assert run_data.AA_processing_dir == first
    assert os.path.isdir(first)


def test_accession_file_run_has_its_processing_dir_at_build_time(spec, tmp_path,
                                                                 listing):
    """The `-a` path was never broken; this pins it so a fix here can't regress it."""
    args = make_args(ncbi_accessions=listing("accs.txt", ["GCF_000008865.2"]),
                     target_pfams_file=listing("t.txt", ["PF90001"]),
                     output_dir=str(tmp_path / "out"))
    out_dir, work_dir = setup_output_dir(args, spec)
    run_data = build_run_data(args, spec, out_dir, work_dir)

    assert os.path.isdir(run_data.ncbi_processing_dir)


def test_reference_data_is_fetched_for_every_ncbi_bearing_input(spec, monkeypatch):
    """
    These subcommands only ever resolved *paths* into the data dirs; nothing fetched
    the assets. A machine that had never run the main program would fail on a missing
    NCBI/GTDB table rather than downloading it.
    """
    import gtotree.utils.ncbi.get_ncbi_assembly_data as ncbi_mod
    import gtotree.utils.gtdb.get_gtdb_data as gtdb_mod

    calls = []
    monkeypatch.setattr(ncbi_mod, "get_ncbi_assembly_data",
                        lambda *a, **kw: calls.append("ncbi"))
    monkeypatch.setattr(gtdb_mod, "get_gtdb_data", lambda *a, **kw: calls.append("gtdb"))

    # a GTDB -w selection still yields NCBI accessions to download, so both are needed
    calls.clear()
    ensure_reference_data(make_args(wanted_ref_tax="Nitrospirota", source="gtdb"), spec)
    assert sorted(calls) == ["gtdb", "ncbi"]

    calls.clear()
    ensure_reference_data(make_args(ncbi_accessions="accs.txt"), spec)
    assert calls == ["ncbi"]

    # local files only: neither table is touched
    calls.clear()
    ensure_reference_data(make_args(amino_acid_files="aa.txt"), spec)
    assert calls == []
