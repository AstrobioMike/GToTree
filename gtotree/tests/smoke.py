import argparse
import json
import shutil
import sys
from contextlib import ExitStack
from importlib import resources
from pathlib import Path

from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from gtotree.cli.parser import parser
from gtotree.main import main as run_gtotree
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

    not_found = run_data.get("num_accs_not_found", 0)
    if not_found:
        report_message(
            f"Smoke test FAILED: {not_found} accessions not found (expected 0)", "red")
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


def _cleanup(cwd, output_dir):
    listing = cwd / LISTING_NAME
    if listing.exists():
        listing.unlink()
    if output_dir.exists():
        shutil.rmtree(output_dir, ignore_errors=True)


def build_parser(parent_subparsers=None):

    desc = ("This program runs an end-to-end smoke test of the installed GToTree "
            "environment against bundled mock amino-acid files and a mock HMM, "
            "verifying a tree is produced. It takes no arguments.")

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
    add_help(optional)
    add_version_arg(optional)

    return parser_


def main():
    # parse (handles -h/-v via the shared help machinery); the smoke test itself
    # takes no arguments, so anything parsed is discarded
    build_parser().parse_args()
    sys.exit(run_smoke_test())


def run_smoke_test(argv=None):
    cwd = Path.cwd()
    listing = cwd / LISTING_NAME
    output_dir = cwd / OUTPUT_NAME

    pkg = resources.files(DATA_PKG)

    with ExitStack() as stack:
        # Resolve packaged data to real filesystem paths. On a conda (unpacked)
        # install as_file() yields the real package path -- so the run info
        # displays the package location as provenance. On a zipped install it
        # yields a temp copy valid only inside this ExitStack, which is why the
        # whole run happens in-block.
        aa_paths = [
            stack.enter_context(resources.as_file(pkg / name))
            for name in AA_FILES
        ]
        hmm_path = stack.enter_context(resources.as_file(pkg / HMM_FILE))

        # Listing lives in cwd; its contents point at the package FASTAs.
        # --- SEAM: this assumes preprocess_genomes READS the amino-acid files
        # in place. If it writes anything beside them, point -A at copies in the
        # output dir instead (materialize the FASTAs into a writable dir first).
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
        except Exception as e:
            report_message(f"Smoke test errored: {e}", "red")
            _cleanup(cwd, output_dir)          # leave nothing behind on error
            return 1

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(run_smoke_test())
