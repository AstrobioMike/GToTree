"""
The one place `gtt search-pfams` and `gtt search-kos` differ.

Everything the two subcommands do differently is a value or a callable on a
`TargetSearchSpec`; everything else is shared code in the sibling modules. The point is
that adding a third target type later should mean writing one spec, not a third copy of
the driver.

The `run_data` attribute names are strings rather than direct references because the
existing search workers, counts-table writers, and hit-combining helpers all read them
off `RunData` by name (`pfam_results_dir` vs `ko_results_dir`, `found_pfam_targets` vs
`found_ko_targets`, and so on). Naming them here lets the shared driver drive the
existing, unmodified helpers for either target type.
"""

from dataclasses import dataclass, field
from typing import Callable, List, Optional


@dataclass(frozen=True)
class TargetSearchSpec:
    """
    Describes one target type end to end.

    Attributes fall into four groups:

      * identity/CLI    -- subcommand name, flags, and the words used in output
      * run_data wiring -- which RunData fields this target type reads and writes
      * behavior        -- the callables that resolve targets, search, and write outputs
      * environment     -- the external dependency and managed dataset it needs
    """

    # --- identity and CLI ------------------------------------------------------
    subcommand: str                  # e.g. "search-pfams"
    target_label: str                # singular, for prose: "Pfam"
    target_label_plural: str         # "Pfams"
    targets_flag: str                # "-p"
    targets_flag_long: str           # "--target-pfams"
    targets_dest: str                # "target_pfams_file"
    default_output_dir: str          # "gtt-pfam-search-output"
    example_target: str              # shown in help, e.g. "PF00789"

    # --- run_data wiring -------------------------------------------------------
    results_dir_attr: str            # "pfam_results_dir"
    results_dir_rel_attr: str        # "pfam_results_dir_rel"
    tmp_results_dir_attr: str        # "tmp_pfam_results_dir"
    targets_file_attr: str           # "target_pfams_file"
    total_targets_attr: str          # "total_pfam_targets"
    found_targets_attr: str          # "found_pfam_targets"
    failed_targets_attr: str         # "failed_pfam_targets"
    searching_done_attr: str         # "additional_pfam_searching_done"
    search_done_flag: str            # GenomeData attr: "pfam_search_done"
    search_failed_flag: str

    # --- output layout ---------------------------------------------------------
    # subdirectories created under the (flattened) output directory
    result_subdirs: List[str]
    # files/dirs whose presence proves the target-collection stage really finished,
    # so a `--resume` doesn't trust a stale done-flag over what's on disk
    target_stage_artifacts: List[str]
    hit_seqs_subdir: str             # "pfam-hit-seqs"
    counts_filename: str             # "pfam-hit-counts.tsv"
    failed_targets_filename: str     # "failed-pfam-targets.txt"

    # --- behavior --------------------------------------------------------------
    # (run_data) -> run_data; resolves the target set and writes the profile assets
    collect_targets: Callable
    # which SearchPlan flag turns this target type's per-genome search on. The fused
    # preprocess-then-search worker in main_stages/processing_genomes.py already
    # dispatches on do_pfam/do_ko, so naming the flag here is enough to drive the
    # existing, unmodified worker for either type
    plan_flag: str
    # (run_data) -> None; writes the hit-counts matrix
    write_counts_table: Callable
    # (target_ids, tmp_results_dir, out_base) -> None
    combine_hits: Callable
    # (force_update=False) -> path; ensures the managed dataset is present
    ensure_data: Callable
    # (run_data) -> str or None; a short version string for the run-info block.
    # Called before run_data exists (the data is fetched up front), so it must
    # accept None and read the version off disk rather than off the run.
    describe_data_version: Optional[Callable] = None

    # --- environment -----------------------------------------------------------
    # external binaries this target type needs on PATH, beyond genome processing
    required_binaries: List[str] = field(default_factory=list)
    # environment variables that must be set for this target type
    required_env_vars: List[str] = field(default_factory=list)
    # True when the search goes through pyhmmer and benefits from a one-time hmmpress
    presses_profiles: bool = False
    # RunData attr holding the combined HMM to press (only when presses_profiles)
    combined_hmm_attr: Optional[str] = None

    def results_dir(self, run_data):
        return getattr(run_data, self.results_dir_attr)

    def tmp_results_dir(self, run_data):
        return getattr(run_data, self.tmp_results_dir_attr)

    def found_targets(self, run_data):
        return getattr(run_data, self.found_targets_attr) or []

    def failed_targets(self, run_data):
        return getattr(run_data, self.failed_targets_attr) or []

    def total_targets(self, run_data):
        return getattr(run_data, self.total_targets_attr) or 0

    def targets_file(self, args):
        return getattr(args, self.targets_dest, None)


def _pfam_spec():
    from gtotree.utils.pfam.pfam_handling import get_additional_pfam_targets
    from gtotree.utils.pfam.get_pfam_data import get_pfam_data, get_stored_pfam_version
    from gtotree.utils.pfam.additional_pfam_searching import (write_pfam_counts_table,
                                                              combine_all_pfam_hits)

    def describe_version(run_data=None):
        import os
        pfam_dir = os.environ.get("Pfam_data_dir")
        if not pfam_dir:
            return None
        return get_stored_pfam_version(pfam_dir)

    return TargetSearchSpec(
        subcommand="search-pfams",
        target_label="Pfam",
        target_label_plural="Pfams",
        targets_flag="-p",
        targets_flag_long="--target-pfams",
        targets_dest="target_pfams_file",
        default_output_dir="gtt-pfam-search-output",
        example_target="PF00789",

        results_dir_attr="pfam_results_dir",
        results_dir_rel_attr="pfam_results_dir_rel",
        tmp_results_dir_attr="tmp_pfam_results_dir",
        targets_file_attr="target_pfams_file",
        total_targets_attr="total_pfam_targets",
        found_targets_attr="found_pfam_targets",
        failed_targets_attr="failed_pfam_targets",
        searching_done_attr="additional_pfam_searching_done",
        search_done_flag="pfam_search_done",
        search_failed_flag="pfam_search_failed",

        result_subdirs=["info", "individual-genome-results", "pfam-hit-seqs",
                        "target-pfam-profiles"],
        target_stage_artifacts=["target-pfam-profiles/all-pfam-targets.hmm"],
        hit_seqs_subdir="pfam-hit-seqs",
        counts_filename="pfam-hit-counts.tsv",
        failed_targets_filename="failed-pfam-targets.txt",

        collect_targets=get_additional_pfam_targets,
        plan_flag="do_pfam",
        write_counts_table=write_pfam_counts_table,
        combine_hits=combine_all_pfam_hits,
        ensure_data=get_pfam_data,
        describe_data_version=describe_version,

        required_binaries=[],
        required_env_vars=["Pfam_data_dir"],
        presses_profiles=True,
        combined_hmm_attr="all_pfam_targets_hmm_path",
    )


def _ko_spec():
    from gtotree.utils.ko.ko_handling import parse_kofamscan_targets
    from gtotree.utils.ko.get_kofamscan_data import get_kofamscan_data
    from gtotree.utils.ko.additional_ko_searching import (write_ko_counts_table,
                                                          combine_all_ko_hits,
                                                          write_out_failed_ko_targets)

    def collect(run_data):
        run_data = parse_kofamscan_targets(run_data)
        write_out_failed_ko_targets(run_data)
        return run_data

    return TargetSearchSpec(
        subcommand="search-kos",
        target_label="KO",
        target_label_plural="KOs",
        targets_flag="-K",
        targets_flag_long="--target-kos",
        targets_dest="target_kos_file",
        default_output_dir="gtt-ko-search-output",
        example_target="K01601",

        results_dir_attr="ko_results_dir",
        results_dir_rel_attr="ko_results_dir_rel",
        tmp_results_dir_attr="tmp_ko_results_dir",
        targets_file_attr="target_kos_file",
        total_targets_attr="total_ko_targets",
        found_targets_attr="found_ko_targets",
        failed_targets_attr="failed_ko_targets",
        searching_done_attr="additional_ko_searching_done",
        search_done_flag="ko_search_done",
        search_failed_flag="ko_search_failed",

        result_subdirs=["individual-genome-results", "ko-hit-seqs",
                        "target-ko-profiles"],
        target_stage_artifacts=["target-kos.tsv", "target-ko-profiles"],
        hit_seqs_subdir="ko-hit-seqs",
        counts_filename="ko-hit-counts.tsv",
        failed_targets_filename="failed-ko-targets.txt",

        collect_targets=collect,
        plan_flag="do_ko",
        write_counts_table=write_ko_counts_table,
        combine_hits=combine_all_ko_hits,
        ensure_data=get_kofamscan_data,
        describe_data_version=None,
        required_binaries=["exec_annotation"],
        required_env_vars=["KO_data_dir"],
        presses_profiles=False,
        combined_hmm_attr=None,
    )


# Built lazily: each builder imports its target type's modules (pandas, pyhmmer,
# the Pfam/KO data handlers), and the `gtt` dispatcher imports every subcommand
# module to construct the parser tree. Building both eagerly at import time would
# make every `gtt` invocation pay for both.
_BUILDERS = {"pfam": _pfam_spec, "ko": _ko_spec}
_CACHE = {}


def get_spec(name):
    """
    Return the spec for 'pfam' or 'ko', building it on first use.
    """
    if name not in _BUILDERS:
        raise KeyError(f"unknown target-search spec: {name!r}")
    if name not in _CACHE:
        _CACHE[name] = _BUILDERS[name]()
    return _CACHE[name]
