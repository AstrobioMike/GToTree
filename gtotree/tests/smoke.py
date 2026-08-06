import argparse
import json
import shutil
import sys
from contextlib import ExitStack
from importlib import resources
from pathlib import Path

from gtotree.cli.common import (CustomRichHelpFormatter, add_help, add_version_arg,
                                wrap_help)
from gtotree.cli.parser import parser
from gtotree.main import main as run_gtotree
from gtotree.utils.context import log_file_var
from gtotree.utils.messaging import report_message, color_text

DATA_PKG = "gtotree.tests.data"
AA_FILES = ["mock-1.faa", "mock-2.faa", "mock-3.faa", "mock-4.faa"]
HMM_FILE = "mock.hmm"
EXPECTED_GENOMES = len(AA_FILES)

LISTING_NAME = "test-amino-acid-files.txt"
OUTPUT_NAME = "test-gtotree-output"


def _verify(output_dir):
    ok = True

    run_data_path = output_dir / "run-files" / "run-data.json"
    if not run_data_path.exists():
        report_message(f"Smoke test FAILED: no run-data.json at {run_data_path}", "red")
        return False

    try:
        run_data = json.loads(run_data_path.read_text())
    except (OSError, json.JSONDecodeError) as e:
        report_message(f"Smoke test FAILED: could not read run-data.json: {e}", "red")
        return False

    genomes = run_data.get("all_input_genomes") or []
    if len(genomes) != EXPECTED_GENOMES:
        report_message(
            f"Smoke test FAILED: expected {EXPECTED_GENOMES} genomes, "
            f"found {len(genomes)}", "red")
        ok = False

    removed = [g for g in genomes if g.get("removed")]
    if removed:
        detail = ", ".join(
            f"{g.get('id')} ({g.get('reason_removed') or 'no reason recorded'})"
            for g in removed)
        report_message(
            f"Smoke test FAILED: {len(removed)} genome(s) were dropped: {detail}", "red")
        ok = False

    tree = OUTPUT_NAME + "/test-gtotree-output.tre"
    tree_path = Path(tree)
    if not tree_path.exists() or tree_path.stat().st_size == 0:
        report_message(f"Smoke test FAILED: no tree produced at {tree_path}", "red")
        ok = False

    if ok:

        print()
        print(color_text(f"{'='*82}", "yellow"))
        report_message(f"{' '*26}GToTree smoke test passed!", "yellow", newline=False)
        print(color_text(f"{'='*82}", "yellow"))
        print()

    return ok


def _cleanup(cwd, output_dir, keep=False):

    listing = cwd / LISTING_NAME
    if listing.exists():
        listing.unlink()

    if keep:
        report_message(f"Test output left in place at:", ii="    ", newline=False, color=None)
        report_message(f"{output_dir}", "yellow", ii="      ", newline=False, trailing_newline=True)
        return

    if output_dir.exists():
        shutil.rmtree(output_dir, ignore_errors=True)


def build_parser(parent_subparsers=None):

    desc = ("This program runs an end-to-end smoke test of the installed GToTree "
            "environment against bundled mock amino-acid files and a mock HMM.")

    if parent_subparsers is not None:
        parser_ = parent_subparsers.add_parser(
            "test",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser_ = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `gtt test`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    optional = parser_.add_argument_group("Optional Parameters")
    optional.add_argument(
        "-k",
        "--keep",
        action="store_true",
        help=wrap_help("Keep the test output directory (\"" + OUTPUT_NAME + "\") "
                       "after the run instead of removing it"),
    )
    add_help(optional)
    add_version_arg(optional)

    return parser_


def main():

    sys.exit(run_smoke_test())


def run_smoke_test(argv=None):

    cli_args = build_parser().parse_args(sys.argv[1:] if argv is None else argv)

    cwd = Path.cwd()
    listing = cwd / LISTING_NAME
    output_dir = cwd / OUTPUT_NAME

    log_file_token = log_file_var.set(log_file_var.get())

    pkg = resources.files(DATA_PKG)

    with ExitStack() as stack:

        aa_paths = [
            stack.enter_context(resources.as_file(pkg / name))
            for name in AA_FILES
        ]
        hmm_path = stack.enter_context(resources.as_file(pkg / HMM_FILE))

        listing.write_text("\n".join(str(p) for p in aa_paths) + "\n")

        if output_dir.exists():
            shutil.rmtree(output_dir, ignore_errors=True)

        gtt_argv = [
            "-A", LISTING_NAME,
            "-H", str(hmm_path),
            "-o", OUTPUT_NAME,
            "-j", "4",
        ]
        args = parser().parse_args(gtt_argv)


        print()
        print(color_text(f"{'='*82}", "yellow"))
        report_message(f"{' '*25}Running GToTree smoke test...", "yellow", newline=False)
        print(color_text(f"{'='*82}", "yellow"))

        try:
            run_gtotree(args)
            ok = _verify(output_dir)
        except SystemExit as e:
            report_message(f"Smoke test exited early (status {e.code})", "red")
            return 1
        except Exception as e:
            report_message(f"Smoke test errored: {e}", "red")
            return 1
        finally:
            _cleanup(cwd, output_dir, keep=cli_args.keep)
            log_file_var.reset(log_file_token)

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(run_smoke_test())
