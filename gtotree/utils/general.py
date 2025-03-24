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
    parallel_used: bool = False
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
class RunData:
    ncbi_accessions: List[str] = field(default_factory=list)
    ncbi_accessions_done: List[str] = field(default_factory=list)
    ncbi_accs_not_found: List[str] = field(default_factory=list)
    ncbi_accs_not_downloaded: List[str] = field(default_factory=list)
    genbank_files: List[str] = field(default_factory=list)
    genbank_files_full_paths: List[str] = field(default_factory=list)
    genbank_files_done: List[str] = field(default_factory=list)
    fasta_files: List[str] = field(default_factory=list)
    fasta_files_full_paths: List[str] = field(default_factory=list)
    fasta_files_done: List[str] = field(default_factory=list)
    amino_acid_files: List[str] = field(default_factory=list)
    amino_acid_files_full_paths: List[str] = field(default_factory=list)
    amino_acid_files_done: List[str] = field(default_factory=list)
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
    genbank_processing_dir_rel: str = ""
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
    def num_ncbi_accessions(self) -> int:
        return len(self.ncbi_accessions)

    @property
    def num_genbank_files(self) -> int:
        return len(self.genbank_files)

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

    def add_done_genbank_file(self, filepath: str):
        # if filepath not in self.genbank_files_done:
        #     self.genbank_files_done.append(filepath)
        print(f"\n\n{filepath}\n\n")
        print(f"\n\n{self.genbank_files}\n\n")
        print(f"\n\n{self.all_remaining_genomes}\n\n")

        if filepath in self.genbank_files:
            candidate = filepath
        else:
            target = os.path.basename(filepath)
            for f in self.genbank_files:
                if os.path.basename(f) == target:
                    candidate = f
                if os.path.basename(f)[:-3] == target:
                    candidate = f
                    break

        if candidate:
            try:
                self.genbank_files.remove(candidate)
            except:
                self.genbank_files.remove(candidate[:-3])
            try:
                self.all_remaining_genomes.remove(candidate)
            except:
                self.all_remaining_genomes.remove(candidate[:-3])
            self.removed_genbank_files.append(candidate)
            self.all_removed_genomes.append(candidate)

        if candidate not in self.genbank_files_done:
            self.genbank_files_done.append(candidate)

    def add_done_fasta_file(self, filepath: str):
        if filepath not in self.fasta_files_done:
            self.fasta_files_done.append(filepath)

    def add_done_amino_acid_file(self, filepath: str):
        if filepath not in self.amino_acid_files_done:
            self.amino_acid_files_done.append(filepath)

    def add_ncbi_acc_not_downloaded(self, accession: str):
        if accession not in self.ncbi_accs_not_downloaded:
            self.ncbi_accs_not_downloaded.append(accession)

    def remove_ncbi_accession(self, accession: str):
        if accession in self.ncbi_accessions:
            self.ncbi_accessions.remove(accession)
            self.all_remaining_genomes.remove(accession)
            self.removed_ncbi_accessions.append(accession)
            self.all_removed_genomes.append(accession)

    def remove_genbank_file(self, filepath: str):
        if filepath in self.genbank_files:
            candidate = filepath
        else:
            target = os.path.basename(filepath)
            for f in self.genbank_files:
                if os.path.basename(f) == target:
                    candidate = f
                    break

        if candidate:
            self.genbank_files.remove(candidate)
            self.all_remaining_genomes.remove(candidate)
            self.removed_genbank_files.append(candidate)
            self.all_removed_genomes.append(candidate)

    def remove_fasta_file(self, filepath: str):
        if filepath in self.fasta_files:
            candidate = filepath
        else:
            target = os.path.basename(filepath)
            for f in self.fasta_files:
                if os.path.basename(f) == target:
                    candidate = f
                    break

        if candidate:
            self.fasta_files.remove(candidate)
            self.all_remaining_genomes.remove(candidate)
            self.removed_fasta_files.append(candidate)
            self.all_removed_genomes.append(candidate)

    def remove_amino_acid_file(self, filepath: str):
        # if filepath in self.amino_acid_files:
        #     self.amino_acid_files.remove(filepath)
        #     self.all_remaining_genomes.remove(filepath)
        #     self.removed_amino_acid_files.append(filepath)
        #     self.all_removed_genomes.append(filepath)

        if filepath in self.amino_acid_files:
            candidate = filepath
        else:
            target = os.path.basename(filepath)
            for f in self.amino_acid_files:
                if os.path.basename(f) == target:
                    candidate = f
                    break

        if candidate:
            self.amino_acid_files.remove(candidate)
            self.all_remaining_genomes.remove(candidate)
            self.removed_amino_acid_files.append(candidate)
            self.all_removed_genomes.append(candidate)

    def replace_in_list(self, list_attr, old_value, new_value):
        # current_list = getattr(self, list_attr, None)
        # if current_list is None:
        #     raise AttributeError(f"'{list_attr}' is not an attribute of RunData.")
        # try:
        #     index = current_list.index(old_value)
        #     current_list[index] = new_value
        # except ValueError:
        #     print(f"Value '{old_value}' not found in '{list_attr}'. No replacement made.")
        # index = self.all_remaining_genomes.index(old_value)
        # self.all_remaining_genomes[index] = new_value
        current_list = getattr(self, list_attr, None)
        if current_list is None:
            raise AttributeError(f"'{list_attr}' is not an attribute of RunData.")

        replaced = False
        try:
            index = current_list.index(old_value)
            current_list[index] = new_value
            replaced = True
        except ValueError:
            old_basename = os.path.basename(old_value)
            for i, item in enumerate(current_list):
                if os.path.basename(item) == old_basename:
                    current_list[i] = new_value
                    replaced = True
                    break
        if not replaced:
            raise ValueError(f"Value '{old_value}' not found in '{list_attr}'. No replacement made.")

        # Now update in all_remaining_genomes.
        replaced = False
        try:
            index = self.all_remaining_genomes.index(old_value)
            self.all_remaining_genomes[index] = new_value
            replaced = True
        except ValueError:
            old_basename = os.path.basename(old_value)
            for i, item in enumerate(self.all_remaining_genomes):
                if os.path.basename(item) == old_basename:
                    self.all_remaining_genomes[i] = new_value
                    replaced = True
                    break
        if not replaced:
            print(f"Value '{old_value}' not found in 'all_remaining_genomes'. No replacement made.")



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
        run_data.ncbi_accessions = entries_list
        run_data.all_input_genomes.extend(entries_list)
        run_data.all_remaining_genomes.extend(entries_list)

    if args.genbank_files:
        with open(args.genbank_files, "r") as f:
            entries_list = f.read().splitlines()
        run_data.genbank_files = entries_list
        run_data.all_input_genomes.extend(entries_list)
        run_data.all_remaining_genomes.extend(entries_list)

    if args.fasta_files:
        with open(args.fasta_files, "r") as f:
            entries_list = f.read().splitlines()
        run_data.fasta_files = entries_list
        run_data.all_input_genomes.extend(entries_list)
        run_data.all_remaining_genomes.extend(entries_list)

    if args.amino_acid_files:
        with open(args.amino_acid_files, "r") as f:
            entries_list = f.read().splitlines()
        run_data.amino_acid_files = entries_list
        run_data.all_input_genomes.extend(entries_list)
        run_data.all_remaining_genomes.extend(entries_list)

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
    cmd += ["--directory", run_data.run_files_dir]
    with open(snakemake_log, "w") as log_file, tqdm(total=tqdm_jobs, desc="    " + description, ncols = 78) as pbar:
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
            report_snakemake_failure(description, snakemake_log)
            report_early_exit()


def run_prodigal(id, run_data, group = ["ncbi", "fasta", "genbank"]):

    if group == "ncbi":
        in_path = f"{run_data.ncbi_downloads_dir}/{id}_genomic.fna"
        out_path = f"{run_data.ncbi_downloads_dir}/{id}_protein.faa"
    elif group == "genbank":
        in_path = f"{id}"
        out_path = f"{run_data.genbank_processing_dir}/{os.basename(id)}_protein.faa"
    elif group == "fasta":
        report_early_exit(f"    Prodigal not yet implemented for fasta files.")
    else:
        report_early_exit(f"    Prodigal not yet implemented for \"{group}\".")

    prodigal_cmd = [
        "prodigal",
        "-c",
        "-q",
        "-i", f"{in_path}",
        "-a", f"{out_path}.tmp",
    ]

    remove_ast_cmd = f"tr -d '*' < {out_path}.tmp > {out_path}"

    try:
        subprocess.run(prodigal_cmd, stdout=subprocess.DEVNULL)
        subprocess.run(remove_ast_cmd, shell=True)
        shutil.rmtree(f"{out_path}.tmp")
        done = True
    except:
        done = False

    return done


def gunzip_if_needed(path):
    if path.endswith(".gz"):
        gunzipped_path = path[:-3]
        with gzip.open(path, 'rb') as f_in, open(gunzipped_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        return gunzipped_path
    else:
        return path
