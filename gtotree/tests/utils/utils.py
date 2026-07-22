import subprocess
import pytest # type: ignore

from gtotree.utils.general import GenomeData


def run_cli(cmd, **kwargs):
    """Run a CLI command, failing the test with captured output on nonzero exit."""
    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if result.returncode != 0:
        pytest.fail(
            f"CLI failed with exit code {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
    return result


class FakeRunData:
    """
    A minimal stand-in for RunData exposing only the surface the NCBI
    assembly-summary parser and downloader touch. Uses REAL GenomeData objects
    (via GenomeData.from_acc) so found/removed bookkeeping matches production;
    only the container/paths are faked, since building a full RunData requires
    many unrelated fields.

    Keep this in sync with parse_ncbi_assembly_summary / preprocessing_genomes
    if their RunData access surface changes.
    """

    def __init__(self, accessions, tmp_dir, nucleotide_mode=False):
        self.tmp_dir = str(tmp_dir)
        self.run_files_dir = str(tmp_dir)
        self.run_files_dir_rel = str(tmp_dir)
        self.ncbi_processing_dir = str(tmp_dir)
        self.nucleotide_mode = nucleotide_mode
        self.ncbi_sub_table_path = None
        self.ncbi_accs = [GenomeData.from_acc(a) for a in accessions]

    def get_input_ncbi_accs(self):
        return [gd.id for gd in self.ncbi_accs]

    def get_ncbi_accs_for_preprocessing(self):
        return [gd for gd in self.ncbi_accs if gd.acc_was_found and not gd.preprocessing_done and not gd.removed]
