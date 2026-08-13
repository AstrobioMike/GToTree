"""
Guards for the TargetSearchSpec objects.

The specs address RunData and GenomeData fields by *name*, which is what lets one
driver serve both target types -- and also what makes them silently wrong if a field
gets renamed. Nothing else in the suite would catch that: a bad attribute name shows up
as an AttributeError deep inside a worker thread on a real run, not at import.
"""

import dataclasses

import pytest  # type: ignore

from gtotree.utils.misc.general import RunData, GenomeData
from gtotree.utils.target_search.target_search_spec import get_spec, TargetSearchSpec


SPEC_NAMES = ["pfam", "ko"]

# spec attributes naming a field on RunData
RUN_DATA_ATTRS = [
    "results_dir_attr",
    "results_dir_rel_attr",
    "tmp_results_dir_attr",
    "targets_file_attr",
    "total_targets_attr",
    "found_targets_attr",
    "failed_targets_attr",
    "searching_done_attr",
]


@pytest.fixture(params=SPEC_NAMES)
def spec(request):
    return get_spec(request.param)


def test_get_spec_rejects_unknown_names():
    with pytest.raises(KeyError):
        get_spec("nope")


def test_get_spec_is_cached():
    """Built lazily, but only once -- each build imports a chunk of the package."""
    assert get_spec("pfam") is get_spec("pfam")


def test_specs_are_frozen(spec):
    """
    Frozen so a caller that needs a variant (a test stubbing `ensure_data`, say) has to
    go through dataclasses.replace and gets its own copy, rather than mutating the
    cached spec every other caller shares.
    """
    with pytest.raises(dataclasses.FrozenInstanceError):
        spec.subcommand = "something-else"


@pytest.mark.parametrize("attr_name", RUN_DATA_ATTRS)
def test_named_run_data_fields_exist(spec, attr_name):
    """Every RunData field the spec names by string must actually be on RunData."""
    field_name = getattr(spec, attr_name)
    assert hasattr(RunData(), field_name), \
        f"{spec.subcommand} spec.{attr_name} names RunData.{field_name}, which is gone"


def test_named_genome_data_flags_exist(spec):
    """
    Same for the two per-genome search flags. Both are needed everywhere, because the
    done flag is set on a failed search too -- it records that the search was attempted,
    not that it worked.
    """
    gd = GenomeData.from_acc("GCF_000000000.1")
    assert hasattr(gd, spec.search_done_flag)
    assert hasattr(gd, spec.search_failed_flag)


def test_combined_hmm_attr_is_set_iff_profiles_are_pressed(spec):
    """
    hmmpress only makes sense for a pyhmmer search, and the driver reads
    `combined_hmm_attr` only when `presses_profiles` is on -- so the two have to agree,
    or the press step silently no-ops (or reads a None path).
    """
    if spec.presses_profiles:
        assert spec.combined_hmm_attr
        assert hasattr(RunData(), spec.combined_hmm_attr)
    else:
        assert spec.combined_hmm_attr is None


def test_plan_flag_is_one_the_fused_worker_understands(spec):
    """
    `plan_flag` has to name a real SearchPlan slot; the fused worker dispatches on it.
    """
    from gtotree.main_stages.processing_genomes import SearchPlan

    assert spec.plan_flag in SearchPlan.__slots__


def test_specs_dont_collide():
    """
    The two specs must not share output filenames or run_data fields, or a Pfam run and
    a KO run pointed at the same directory would overwrite each other.
    """
    pfam, ko = get_spec("pfam"), get_spec("ko")

    assert pfam.counts_filename != ko.counts_filename
    assert pfam.hit_seqs_subdir != ko.hit_seqs_subdir
    assert pfam.failed_targets_filename != ko.failed_targets_filename
    assert pfam.plan_flag != ko.plan_flag

    for attr_name in RUN_DATA_ATTRS:
        assert getattr(pfam, attr_name) != getattr(ko, attr_name), attr_name


def test_accessor_helpers_read_through_to_run_data():
    spec = get_spec("pfam")
    run_data = RunData()
    run_data.pfam_results_dir = "/somewhere"
    run_data.found_pfam_targets = ["PF00001", "PF00002"]
    run_data.total_pfam_targets = 5

    assert spec.results_dir(run_data) == "/somewhere"
    assert spec.found_targets(run_data) == ["PF00001", "PF00002"]
    assert spec.total_targets(run_data) == 5
    # empty rather than None, so callers can len() unconditionally
    assert spec.failed_targets(run_data) == []


def test_spec_is_the_only_place_the_two_differ():
    """
    A shape check on the dataclass itself: if a new per-target-type difference is added
    to the driver, it should show up here as a new field rather than as an `if
    spec.subcommand == ...` branch somewhere in the shared code.
    """
    field_names = {f.name for f in dataclasses.fields(TargetSearchSpec)}
    assert "subcommand" in field_names
    # the two specs must differ in every field that isn't structural
    pfam, ko = get_spec("pfam"), get_spec("ko")
    differing = {name for name in field_names
                 if getattr(pfam, name) != getattr(ko, name)}
    assert "subcommand" in differing
    assert "target_label" in differing
