# import inspect
import os
import sys
import gzip
import shutil
import json #type: ignore
import argparse
from dataclasses import dataclass, field, asdict
from typing import List
from tqdm import tqdm # type: ignore
import subprocess
import urllib.request
from pkg_resources import resource_filename
from gtotree.utils.messaging import report_early_exit, report_snakemake_failure


@dataclass
class ToolsUsed:
    prodigal_used: bool = False
    taxonkit_used: bool = False
    gtdb_used: bool = False
    fasttree_used: bool = False
    veryfasttree_used: bool = False
    iqtree_used: bool = False
    universal_SCGs_used: bool = False
    pfam_db_used: bool = False
    kofamscan_used: bool = False


def download_with_tqdm(url, target, filename=None, urlopen=False):
    with tqdm(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc=target, ncols = 90) as t:
        def reporthook(block_num, block_size, total_size):
            if total_size > 0:
                t.total = total_size
            t.update(block_size)
        if not urlopen:
            urllib.request.urlretrieve(url, filename, reporthook=reporthook)
            sys.stdout.write("")
        else:
            dl = urllib.request.urlopen(url, reporthook=reporthook)
            sys.stdout.write("")
            return dl


@dataclass
class GenomeData:
    id: str
    full_path: str
    provided_path: str
    basename: str
    done: bool = False
    final_AA_path: str = ""
    removed: bool = False
    prodigal_used: bool = False
    was_gzipped: bool = False

    @classmethod
    def from_path(cls, path: str):
        full_path = os.path.abspath(path)
        provided_path = path
        basename = os.path.basename(full_path)

        extensions_to_remove = [".gb", ".gbff", ".fasta", ".fna", ".fa", ".faa"]

        id = basename
        if id.lower().endswith(".gz"):
            id = id[:-3]

        for ext in extensions_to_remove:
            if id.lower().endswith(ext):
                id = id[:-len(ext)]
                break

        return cls(id, full_path, provided_path, basename)

    @classmethod
    def from_acc(cls, acc: str):
        full_path = None
        provided_path = None
        basename = None
        id = acc

        return cls(id, full_path, provided_path, basename)

    def mark_done(self, value=True):
        self.done = value

    def mark_removed(self, value=True):
        self.removed = value
        self.final_AA_path = None

    def mark_prodigal_used(self, value=True):
        self.prodigal_used = value

    def mark_was_gzipped(self, value=True):
        self.was_gzipped = value


@dataclass
class RunData:
    input_ncbi_accessions: List[str] = field(default_factory=list)
    remaining_ncbi_accessions: List[str] = field(default_factory=list)
    ncbi_accessions_done: List[str] = field(default_factory=list)
    ncbi_accs_not_done: List[str] = field(default_factory=list)
    ncbi_accs_not_found: List[str] = field(default_factory=list)
    ncbi_accs_not_downloaded: List[str] = field(default_factory=list)
    genbank_files: List[GenomeData] = field(default_factory=list)
    fasta_files: List[GenomeData] = field(default_factory=list)
    amino_acid_files: List[GenomeData] = field(default_factory=list)
    all_input_genomes: List[str] = field(default_factory=list)
    all_remaining_genomes: List[str] = field(default_factory=list)

    removed_ncbi_accessions: List[str] = field(default_factory=list)
    removed_genbank_files: List[str] = field(default_factory=list)
    removed_fasta_files: List[str] = field(default_factory=list)
    removed_amino_acid_files: List[str] = field(default_factory=list)
    all_removed_genomes: List[str] = field(default_factory=list)

    ncbi_sub_table_path: str = ""
    ncbi_downloads_dir: str = ""
    ncbi_downloads_dir_rel: str = ""
    genbank_processing_dir: str = ""
    fasta_processing_dir: str = ""
    AA_processing_dir: str = ""
    run_files_dir: str = ""
    run_files_dir_rel: str = ""
    run_data_path: str = ""
    all_input_genome_AA_files_dir: str = ""
    tmp_dir: str = ""
    log_file: str = ""
    snakemake_logs_dir: str = ""
    snakemake_logs_dir_rel: str = ""

    tools_used: ToolsUsed = field(default_factory=ToolsUsed)

    @property
    def num_input_ncbi_accessions(self) -> int:
        return len(self.input_ncbi_accessions)

    @property
    def num_remaining_ncbi_accessions(self) -> int:
        return len(self.remaining_ncbi_accessions)

    @property
    def num_ncbi_accs_not_done(self) -> int:
        return len(self.ncbi_accs_not_done)

    @property
    def num_genbank_files(self) -> int:
        return len(self.genbank_files)

    @property
    def num_incomplete_genbank_files(self) -> int:
        return len([gf for gf in self.genbank_files if not gf.done and not gf.removed])

    @property
    def num_incomplete_fasta_files(self) -> int:
        return len([gf for gf in self.fasta_files if not gf.done and not gf.removed])

    @property
    def num_incomplete_amino_acid_files(self) -> int:
        return len([gf for gf in self.amino_acid_files if not gf.done and not gf.removed])

    @property
    def num_fasta_files(self) -> int:
        return len(self.fasta_files)

    @property
    def num_amino_acid_files(self) -> int:
        return len(self.amino_acid_files)

    @property
    def num_input_genomes(self) -> int:
        return len(self.all_input_genomes)

    @property
    def num_remaining_genomes(self) -> int:
        return len(self.all_remaining_genomes)

    @property
    def num_removed_genomes(self) -> int:
        return len(self.all_removed_genomes)

    def add_done_ncbi_accession(self, accession: str):
        if accession not in self.ncbi_accessions_done:
            self.ncbi_accessions_done.append(accession)
        if accession in self.removed_ncbi_accessions:
            self.removed_ncbi_accessions.remove(accession)
        if accession in self.all_removed_genomes:
            self.all_removed_genomes.remove(accession)
        if accession not in self.all_remaining_genomes:
            self.all_remaining_genomes.append(accession)

    def add_ncbi_acc_not_downloaded(self, accession: str):
        if accession not in self.ncbi_accs_not_downloaded:
            self.ncbi_accs_not_downloaded.append(accession)
        if accession not in self.ncbi_accs_not_done:
            self.ncbi_accs_not_done.append(accession)

    def remove_ncbi_accession(self, accession: str):
        if accession in self.remaining_ncbi_accessions:
            self.remaining_ncbi_accessions.remove(accession)
            self.all_remaining_genomes.remove(accession)
            self.removed_ncbi_accessions.append(accession)
            self.all_removed_genomes.append(accession)

    def add_genbank_file(self, filepath):
        gf = GenomeData.from_path(filepath)
        if gf not in self.genbank_files:
            self.genbank_files.append(gf)

    def add_fasta_file(self, filepath):
        gf = GenomeData.from_path(filepath)
        if gf not in self.fasta_files:
            self.fasta_files.append(gf)

    def add_amino_acid_file(self, filepath):
        gf = GenomeData.from_path(filepath)
        if gf not in self.amino_acid_files:
            self.amino_acid_files.append(gf)

    def incomplete_genbank_files(self) -> List[GenomeData]:
        return [gf.provided_path for gf in self.genbank_files if not gf.done and not gf.removed]

    def incomplete_fasta_files(self) -> List[GenomeData]:
        return [gf.provided_path for gf in self.fasta_files if not gf.done and not gf.removed]

    def incomplete_amino_acid_files(self) -> List[GenomeData]:
        return [gf.provided_path for gf in self.amino_acid_files if not gf.done and not gf.removed]

    def failed_genbank_files(self) -> List[GenomeData]:
        return [gf.provided_path for gf in self.genbank_files if not gf.done and gf.removed]

    def failed_fasta_files(self) -> List[GenomeData]:
        return [gf.provided_path for gf in self.fasta_files if not gf.done and gf.removed]

    def failed_amino_acid_files(self) -> List[GenomeData]:
        return [gf.provided_path for gf in self.amino_acid_files if not gf.done and gf.removed]

    def any_incomplete_genbank_files(self) -> bool:
        return any(not gf.done and not gf.removed for gf in self.genbank_files)

    def any_incomplete_fasta_files(self) -> bool:
        return any(not gf.done and not gf.removed for gf in self.fasta_files)

    def any_incomplete_amino_acid_files(self) -> bool:
        return any(not gf.done and not gf.removed for gf in self.amino_acid_files)

    def genbank_files_with_prodigal_used(self) -> List[GenomeData]:
        return [gf for gf in self.genbank_files if gf.prodigal_used]

    def add_done_fasta_file(self, filepath: str):
        if filepath not in self.fasta_files_done:
            self.fasta_files_done.append(filepath)

    def add_done_amino_acid_file(self, filepath: str):
        if filepath not in self.amino_acid_files_done:
            self.amino_acid_files_done.append(filepath)

    def remove_genbank_file(self, filepath: str):
        candidate = None
        full_path = os.path.abspath(filepath)

        for gf in self.genbank_files:
            if gf.full_path == full_path:
                candidate = gf.provided_path
                break

        if candidate is None:
            target = os.path.basename(filepath)
            for gf in self.genbank_files:
                if gf.basename == target:
                    candidate = gf.provided_path
                    break

        if candidate:
            if candidate in self.all_remaining_genomes:
                self.all_remaining_genomes.remove(candidate)
            if full_path in self.all_remaining_genomes:
                self.all_remaining_genomes.remove(full_path)
            if candidate not in self.all_removed_genomes:
                self.all_removed_genomes.append(candidate)
            if candidate not in self.removed_genbank_files:
                self.removed_genbank_files.append(candidate)


    def remove_fasta_file(self, filepath: str):
        candidate = None
        full_path = os.path.abspath(filepath)

        for gf in self.fasta_files:
            if gf.full_path == full_path:
                candidate = gf.provided_path
                break

        if candidate is None:
            target = os.path.basename(filepath)
            for gf in self.fasta_files:
                if gf.basename == target:
                    candidate = gf.provided_path
                    break

        if candidate:
            if candidate in self.all_remaining_genomes:
                self.all_remaining_genomes.remove(candidate)
            if full_path in self.all_remaining_genomes:
                self.all_remaining_genomes.remove(full_path)
            if candidate not in self.all_removed_genomes:
                self.all_removed_genomes.append(candidate)
            if candidate not in self.removed_fasta_files:
                self.removed_fasta_files.append(candidate)


    def remove_amino_acid_file(self, filepath: str):
        candidate = None
        full_path = os.path.abspath(filepath)

        for gf in self.amino_acid_files:
            if gf.full_path == full_path:
                candidate = gf.provided_path
                break

        if candidate is None:
            target = os.path.basename(filepath)
            for gf in self.amino_acid_files:
                if gf.basename == target:
                    candidate = gf.provided_path
                    break

        if candidate:
            if candidate in self.all_remaining_genomes:
                self.all_remaining_genomes.remove(candidate)
            if full_path in self.all_remaining_genomes:
                self.all_remaining_genomes.remove(full_path)
            if candidate not in self.all_removed_genomes:
                self.all_removed_genomes.append(candidate)
            if candidate not in self.removed_amino_acid_files:
                self.removed_amino_acid_files.append(candidate)



    # def __setattr__(self, name, value):
    #     """Debug when an attribute is changed."""
    #     stack = inspect.stack()
    #     caller = stack[1]  # Get the function that triggered the change
    #     print(f"Setting `{name}` to `{value}` in {caller.filename}:{caller.lineno}, function {caller.function}")
    #     super().__setattr__(name, value)


def populate_run_data(args):
    run_data = RunData()

    if args.ncbi_accessions:
        with open(args.ncbi_accessions, "r") as f:
            entries_list = f.read().splitlines()
        run_data.input_ncbi_accessions = entries_list
        run_data.remaining_ncbi_accessions = entries_list
        run_data.all_input_genomes.extend(entries_list)
        run_data.all_remaining_genomes.extend(entries_list)

    if args.genbank_files:
        with open(args.genbank_files, "r") as f:
            entries_list = f.read().splitlines()
        run_data.genbank_files = [GenomeData.from_path(entry) for entry in entries_list]
        run_data.all_input_genomes.extend([gf.full_path for gf in run_data.genbank_files])
        run_data.all_remaining_genomes.extend([gf.full_path for gf in run_data.genbank_files])

    if args.fasta_files:
        with open(args.fasta_files, "r") as f:
            entries_list = f.read().splitlines()
        run_data.fasta_files = [GenomeData.from_path(entry) for entry in entries_list]
        for gf in run_data.fasta_files:
            gf.prodigal_used = True
        run_data.all_input_genomes.extend([gf.full_path for gf in run_data.fasta_files])
        run_data.all_remaining_genomes.extend([gf.full_path for gf in run_data.fasta_files])

    if args.amino_acid_files:
        with open(args.amino_acid_files, "r") as f:
            entries_list = f.read().splitlines()
        run_data.amino_acid_files = [GenomeData.from_path(entry) for entry in entries_list]
        run_data.all_input_genomes.extend([gf.full_path for gf in run_data.amino_acid_files])
        run_data.all_remaining_genomes.extend([gf.full_path for gf in run_data.amino_acid_files])

    run_data.run_files_dir = args.run_files_dir
    run_data.run_files_dir_rel = args.run_files_dir_rel
    run_data.run_data_path = run_data.run_files_dir + "/genome-data.json"

    return run_data


def write_args(args):
    args_path = args.tmp_dir + "/args.json"
    with open(args_path, "w") as f:
        json.dump(args.__dict__, f)
    return args_path


def read_args(args_path):
    with open(args_path, "r") as f:
        args_dict = json.load(f)
    args = argparse.Namespace(**args_dict)
    return args


def write_run_data(run_data):
    with open(run_data.run_data_path, "w") as f:
        json.dump(asdict(run_data), f, indent=2)


def read_run_data(path):
    try:
        with open(path, "r") as f:
            run_data_dict = json.load(f)

        if "genbank_files" in run_data_dict:
            run_data_dict["genbank_files"] = [GenomeData(**gf) if isinstance(gf, dict) else gf for gf in run_data_dict["genbank_files"]]
        if "fasta_files" in run_data_dict:
            run_data_dict["fasta_files"] = [GenomeData(**gf) if isinstance(gf, dict) else gf for gf in run_data_dict["fasta_files"]]
        if "amino_acid_files" in run_data_dict:
            run_data_dict["amino_acid_files"] = [GenomeData(**gf) if isinstance(gf, dict) else gf for gf in run_data_dict["amino_acid_files"]]

        if "tools_used" in run_data_dict and run_data_dict["tools_used"] is not None:
            run_data_dict["tools_used"] = ToolsUsed(**run_data_dict["tools_used"])
        else:
            run_data_dict["tools_used"] = ToolsUsed()

        run_data = RunData(**run_data_dict)
        return run_data

    except FileNotFoundError:
        pass


def get_snakefile_path(basename):
    return resource_filename("gtotree", f"smk/{basename}")


def touch(path):
    with open(path, 'a'):
        os.utime(path, None)


def run_snakemake(cmd, tqdm_jobs, run_data, description, print_lines=False):
    print("")
    num_finished = 0 # counting this so it doesn't flip the progress bar when it counts the final "all" rule
    snakemake_log = f"{run_data.snakemake_logs_dir}/{description.replace(' ', '-').lower()}.log"
    snakemake_log_dir_rel = f"{run_data.snakemake_logs_dir_rel}/{description.replace(' ', '-').lower()}.log"
    cmd += ["--directory", run_data.run_files_dir]
    bar_format = "      {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]"
    with open(snakemake_log, "w") as log_file, tqdm(total = tqdm_jobs, bar_format = bar_format, ncols = 76) as pbar:
    # with open(snakemake_log, "w") as log_file, tqdm(total=tqdm_jobs, desc="      " + description, ncols = 74) as pbar:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )
        for line in iter(process.stdout.readline, ""):
            log_file.write(line)
            log_file.flush()

            if "Finished job" in line:
                num_finished += 1
                if num_finished <= tqdm_jobs:
                    pbar.update(1)
            if print_lines:
                print(line, end="")

        process.wait()

        if process.returncode != 0:
            report_snakemake_failure(description, snakemake_log_dir_rel)
            report_early_exit()


def run_prodigal(id, run_data, full_inpath = None, group = ["ncbi", "fasta", "genbank"]):

    if group == "ncbi":
        in_path = f"{run_data.ncbi_downloads_dir}/{id}_genomic.fna"
        out_path = f"{run_data.ncbi_downloads_dir}/{id}_protein.faa"
    elif group == "genbank":
        in_path = f"{run_data.genbank_processing_dir}/{id}.fasta"
        out_path = f"{run_data.genbank_processing_dir}/{id}_protein.faa"
        print(f"\n\n    {out_path}\n\n")
    elif group == "fasta":
        in_path = full_inpath
        out_path = f"{run_data.fasta_processing_dir}/{id}_protein.faa"
    else:
        report_early_exit(f"    Prodigal not yet implemented for \"{group}\".")

    prodigal_cmd = [
        "prodigal",
        "-c",
        "-q",
        "-i", f"{in_path}",
        "-a", f"{out_path}.tmp",
    ]
    print(f"\n\n    {prodigal_cmd}\n\n")

    remove_ast_cmd = f"tr -d '*' < {out_path}.tmp > {out_path}"

    try:
        subprocess.run(prodigal_cmd, stdout=subprocess.DEVNULL)
        subprocess.run(remove_ast_cmd, shell=True)
        os.remove(f"{out_path}.tmp")
        done = True
    except:
        done = False

    if os.path.getsize(out_path) == 0:
        os.remove(out_path)
        done = False

    return done


def gunzip_if_needed(path):
    if path.endswith(".gz"):
        gunzipped_path = path[:-3]
        with gzip.open(path, 'rb') as f_in, open(gunzipped_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        return gunzipped_path, True
    else:
        return path, False


def remove_file_if_exists(path):
    try:
        os.remove(path)
    except FileNotFoundError:
        pass
