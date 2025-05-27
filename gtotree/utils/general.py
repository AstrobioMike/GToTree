# import inspect
import os
import sys
import gzip
import shutil
import json
import argparse
import pandas as pd
from dataclasses import dataclass, field, asdict
from typing import List
from tqdm import tqdm # type: ignore
import subprocess
import urllib.request
from datetime import datetime
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
    pfam_db_used: bool = False
    kofamscan_used: bool = False
    universal_SCGs_used: bool = False


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


def download_and_gunzip(url, target):

    try:
        with urllib.request.urlopen(url) as resp, gzip.GzipFile(fileobj=resp) as gzipped, open(target, "wb") as out_f:
            shutil.copyfileobj(gzipped, out_f)
        return True

    except:
        return False


@dataclass
class GenomeData:
    id: str
    source: str
    full_path: str
    provided_path: str
    basename: str
    taxid: str = None
    preprocessing_done: bool = False
    preprocessing_failed: bool = False
    final_AA_path: str = ""
    prodigal_used: bool = False
    was_gzipped: bool = False
    acc_was_found: bool = None
    acc_was_downloaded: bool = None
    mapping: str = None
    hmm_search_done: bool = False
    hmm_search_failed: bool = None
    extract_seqs_failed: bool = None
    ko_search_done: bool = False
    ko_search_failed: bool = False
    pfam_search_done: bool = False
    pfam_search_failed: bool = False
    num_genes: int = None
    num_SCG_hits: int = None
    num_unique_SCG_hits: int = None
    num_SCG_hits_after_filtering: int = None
    reason_removed: str = None
    removed: bool = False

    @classmethod
    def from_path(cls, path: str, source: str):
        full_path = os.path.abspath(path)
        provided_path = path
        basename = os.path.basename(full_path)
        source = source

        extensions_to_remove = [".gb", ".gbff", ".fasta", ".fna", ".fa", ".faa"]

        id = basename
        if id.lower().endswith(".gz"):
            id = id[:-3]

        for ext in extensions_to_remove:
            if id.lower().endswith(ext):
                id = id[:-len(ext)]
                break

        return cls(id=id,
                   source=source,
                   full_path=full_path,
                   provided_path=provided_path,
                   basename=basename)

    @classmethod
    def from_acc(cls, acc: str):
        full_path = None
        provided_path = None
        basename = acc
        id = acc
        source = "accession"

        return cls(id=id,
                   source=source,
                   full_path=full_path,
                   provided_path=provided_path,
                   basename=basename)

    def mark_preprocessing_done(self, value=True):
        self.preprocessing_done = value

    def mark_removed(self, reason, value=True):
        self.removed = value
        self.reason_removed = reason
        self.final_AA_path = None

    def mark_prodigal_used(self, value=True):
        self.prodigal_used = value

    def mark_was_gzipped(self, value=True):
        self.was_gzipped = value

    def mark_hmm_search_done(self, value=True):
        self.hmm_search_done = value
        self.hmm_search_failed = False

    def mark_hmm_search_failed(self, value=True):
        self.hmm_search_failed = value

    def mark_extract_seqs_failed(self, value=True):
        self.extract_seqs_failed = value

    def mark_ko_search_done(self, value=True):
        self.ko_search_done = value

    def mark_ko_search_failed(self, value=True):
        self.ko_search_done = True
        self.ko_search_failed = value

    def mark_pfam_search_done(self, value=True):
        self.pfam_search_done = value

    def mark_pfam_search_failed(self, value=True):
        self.pfam_search_done = True
        self.pfam_search_failed = value


@dataclass
class SCGset:
    id: str
    remaining: bool = None
    gene_length_filtered: bool = None
    aligned: bool = None
    trimmed: bool = None
    ready_for_cat: bool = None
    reason_removed: str = None
    removed: bool = None

    @classmethod
    def from_id(cls, id: str):
        return cls(id, remaining=True, gene_length_filtered=False, removed=False)

    def mark_removed(self, reason, value=True):
        self.removed = value
        self.reason_removed = reason
        self.remaining = False

    def mark_gene_length_filtered(self, value=True):
        self.gene_length_filtered = value


@dataclass
class RunData:
    ncbi_accs: List[GenomeData] = field(default_factory=list)
    genbank_files: List[GenomeData] = field(default_factory=list)
    fasta_files: List[GenomeData] = field(default_factory=list)
    amino_acid_files: List[GenomeData] = field(default_factory=list)
    all_input_genomes: List[GenomeData] = field(default_factory=list)
    SCG_targets: List[SCGset] = field(default_factory=list)

    start_time: str = None
    ncbi_sub_table_path: str = ""
    ncbi_downloads_dir: str = ""
    ncbi_downloads_dir_rel: str = ""
    genbank_processing_dir: str = ""
    fasta_processing_dir: str = ""
    AA_processing_dir: str = ""
    output_dir: str = ""
    output_dir_rel: str = ""
    run_files_dir: str = ""
    run_files_dir_rel: str = ""
    individual_gene_alignments_dir: str = ""
    individual_gene_alignments_dir_rel: str = ""
    run_data_path: str = ""
    hmm_path: str = ""
    mapping_file_path: str = ""
    mapping_dict: dict = field(default_factory=dict)
    initial_mapping_IDs_from_user: List[str] = field(default_factory=list)
    ready_genome_AA_files_dir: str = ""
    hmm_results_dir: str = ""
    found_SCG_seqs_dir: str = ""
    num_SCG_targets_removed: int = 0
    tmp_dir: str = ""
    log_file: str = ""
    logs_dir: str = ""
    logs_dir_rel: str = ""
    gtotree_logs_dir: str = ""
    num_hmm_cpus: str = ""
    best_hit_mode: bool = False
    seq_length_cutoff: float = None
    SCG_hits_filtered: bool = False
    genomes_filtered_for_min_SCG_hits: bool = False
    all_SCG_sets_aligned: bool = False
    updating_headers: bool = False
    headers_updated: bool = False
    use_muscle_super5: bool = False
    num_muscle_threads: int = 5
    nucleotide_mode: bool = False
    concatenated_alignment_path: str = ""
    final_alignment_path: str = ""
    final_alignment_length: int = 0
    original_tree_path: str = ""
    final_tree_path: str = ""
    target_kos_file: str = None
    total_ko_targets: int = 0
    target_kos_tsv: str = None
    target_ko_profiles_dir: str = None
    target_pfams_file: str = None
    total_pfam_targets: int = 0
    additional_pfam_searching_done: bool = False
    additional_ko_searching_done: bool = False

    pfam_dict: dict = field(default_factory=dict)
    pfam_results_dir: str = ""
    pfam_results_dir_rel: str = ""
    wanted_pfam_targets: List[str] = field(default_factory=list)
    found_pfam_targets: List[str] = field(default_factory=list)
    failed_pfam_targets: List[str] = field(default_factory=list)

    ko_results_dir: str = ""
    ko_results_dir_rel: str = ""
    tmp_ko_results_dir: str = ""
    wanted_ko_targets: List[str] = field(default_factory=list)
    found_ko_targets: List[str] = field(default_factory=list)
    failed_ko_targets: List[str] = field(default_factory=list)

    tools_used: ToolsUsed = field(default_factory=ToolsUsed)

    @property
    def num_incomplete_genbank_files(self) -> int:
        return len([gd for gd in self.genbank_files if not gd.preprocessing_done and not gd.removed])

    @property
    def num_incomplete_fasta_files(self) -> int:
        return len([gd for gd in self.fasta_files if not gd.preprocessing_done and not gd.removed])

    @property
    def num_incomplete_amino_acid_files(self) -> int:
        return len([gd for gd in self.amino_acid_files if not gd.preprocessing_done and not gd.removed])

    def num_accs_not_found(self) -> int:
        return len([gd for gd in self.ncbi_accs if gd.acc_not_found is True])

    def get_all_SCG_targets(self) -> List[SCGset]:
        return [scg for scg in self.SCG_targets]

    def get_all_SCG_targets_remaining(self) -> List[SCGset]:
        return [scg for scg in self.SCG_targets if scg.remaining and not scg.removed]

    def get_all_removed_SCG_targets(self) -> List:
        return [scg.id for scg in self.SCG_targets if scg.removed]

    def get_all_SCG_targets_remaining_but_not_filtered(self) -> List[SCGset]:
        return [scg for scg in self.SCG_targets if scg.remaining and not scg.gene_length_filtered]

    def get_all_SCG_targets_remaining_but_not_aligned(self) -> List[SCGset]:
        return [scg for scg in self.SCG_targets if scg.remaining and not scg.aligned]

    def get_all_SCG_targets_ready_for_concatenation(self) -> List[SCGset]:
        return [scg for scg in self.SCG_targets if scg.ready_for_cat and not scg.removed]

    def update_all_input_genomes(self):
        self.all_input_genomes = []
        self.all_input_genomes.extend(self.ncbi_accs)
        self.all_input_genomes.extend(self.genbank_files)
        self.all_input_genomes.extend(self.fasta_files)
        self.all_input_genomes.extend(self.amino_acid_files)

    def get_all_input_genome_ids(self) -> List[str]:
        return [gd.id for gd in self.all_input_genomes]

    def get_all_preprocessed_genomes(self) -> List[GenomeData]:
        return [gd for gd in self.all_input_genomes if gd.preprocessing_done]

    def get_all_input_genomes_for_hmm_search(self) -> List[GenomeData]:
        return [gd for gd in self.all_input_genomes if gd.preprocessing_done and not gd.hmm_search_done and not gd.removed]

    def get_all_input_genomes_for_ko_search(self) -> List[GenomeData]:
        return [gd for gd in self.all_input_genomes if gd.preprocessing_done and not gd.ko_search_done and not gd.removed]

    def get_all_input_genomes_for_pfam_search(self) -> List[GenomeData]:
        return [gd for gd in self.all_input_genomes if gd.preprocessing_done and not gd.pfam_search_done and not gd.removed]

    def get_all_input_genomes_for_filtering(self) -> List[GenomeData]:
        return [gd for gd in self.all_input_genomes if gd.preprocessing_done and gd.hmm_search_done and not gd.removed]

    def get_all_input_genomes_due_for_SCG_min_hit_filtering(self) -> List[GenomeData]:
        return [gd for gd in self.all_input_genomes if gd.reason_removed == "too few SCG hits" or gd.reason_removed == "too few unique SCG hits"]

    def get_all_input_genome_basenames(self) -> List[str]:
        return [gd.basename for gd in self.all_input_genomes]

    def get_all_input_genome_provided_paths(self) -> List[str]:
        return [gd.provided_path for gd in self.all_input_genomes]

    def get_input_ncbi_accs(self) -> List[str]:
        return [gd.id for gd in self.ncbi_accs]

    def get_input_genbank_ids(self) -> List[str]:
        return [gd.id for gd in self.genbank_files]

    def get_input_fasta_ids(self) -> List[str]:
        return [gd.id for gd in self.fasta_files]

    def get_input_amino_acid_ids(self) -> List[str]:
        return [gd.id for gd in self.amino_acid_files]

    def found_ncbi_accs(self) -> List[GenomeData]:
        return [gd for gd in self.ncbi_accs if gd.acc_was_found]

    def get_ncbi_accs_for_snakemake_preprocessing(self) -> List[GenomeData]:
        return [gd for gd in self.ncbi_accs if gd.acc_was_found and not gd.preprocessing_done and not gd.removed]

    def remaining_ncbi_accs(self) -> List[GenomeData]:
        return [gd.id for gd in self.ncbi_accs if not gd.removed]

    def get_ncbi_accs_not_downloaded(self) -> List[str]:
        return [gd.id for gd in self.ncbi_accs if gd.acc_was_downloaded is False]

    def get_ncbi_accs_not_found(self) -> List[str]:
        return [gd.id for gd in self.ncbi_accs if gd.acc_was_found is False]

    def get_removed_ncbi_accs(self) -> List[str]:
        return [gd.id for gd in self.ncbi_accs if gd.removed]

    def get_all_removed_input_genomes(self) -> List[str]:
        return [gd.id for gd in self.all_input_genomes if gd.removed]

    def get_all_remaining_input_genomes(self) -> List[GenomeData]:
        return [gd for gd in self.all_input_genomes if not gd.removed]

    def get_all_remaining_input_genome_ids(self) -> List[str]:
        return [gd.id for gd in self.all_input_genomes if not gd.removed]

    def get_done_ncbi_accs(self) -> List[GenomeData]:
        return [gd for gd in self.ncbi_accs if gd.preprocessing_done]

    def get_failed_genbank_ids(self) -> List[str]:
        return [gd.id for gd in self.genbank_files if gd.removed]

    def get_failed_genbank_paths(self) -> List[str]:
        return [gd.provided_path for gd in self.genbank_files if gd.removed]

    def get_failed_fasta_ids(self) -> List[str]:
        return [gd.id for gd in self.fasta_files if gd.removed]

    def get_failed_fasta_paths(self) -> List[str]:
        return [gd.provided_path for gd in self.fasta_files if gd.removed]

    def get_failed_amino_acid_ids(self) -> List[str]:
        return [gd.id for gd in self.amino_acid_files if gd.removed and gd.preprocessing_failed]

    def get_failed_amino_acid_paths(self) -> List[str]:
        return [gd.provided_path for gd in self.amino_acid_files if gd.removed and gd.preprocessing_failed]

    def get_failed_hmm_search_paths(self) -> List[str]:
        return [gd.provided_path for gd in self.all_input_genomes if gd.preprocessing_done and not gd.hmm_search_done]

    def get_prodigal_used_genbank_ids(self) -> List[str]:
        return [gd.id for gd in self.genbank_files if gd.prodigal_used]

    def incomplete_genbank_files(self) -> List[GenomeData]:
        return [gd.provided_path for gd in self.genbank_files if not gd.preprocessing_done and not gd.removed]

    def incomplete_fasta_files(self) -> List[GenomeData]:
        return [gd.provided_path for gd in self.fasta_files if not gd.preprocessing_done and not gd.removed]

    def incomplete_amino_acid_files(self) -> List[GenomeData]:
        return [gd.provided_path for gd in self.amino_acid_files if not gd.preprocessing_done and not gd.removed]

    def any_incomplete_genbank_files(self) -> bool:
        return any(not gd.preprocessing_done and not gd.removed for gd in self.genbank_files)

    def any_incomplete_fasta_files(self) -> bool:
        return any(not gd.preprocessing_done and not gd.removed for gd in self.fasta_files)

    def any_incomplete_amino_acid_files(self) -> bool:
        return any(not gd.preprocessing_done and not gd.removed for gd in self.amino_acid_files)

    def genbank_files_with_prodigal_used(self) -> List[GenomeData]:
        return [gd for gd in self.genbank_files if gd.prodigal_used]


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
        run_data.ncbi_accs = [GenomeData.from_acc(entry) for entry in entries_list]

    if args.genbank_files:
        with open(args.genbank_files, "r") as f:
            entries_list = f.read().splitlines()
        run_data.genbank_files = [GenomeData.from_path(entry, "genbank-file") for entry in entries_list]

    if args.fasta_files:
        with open(args.fasta_files, "r") as f:
            entries_list = f.read().splitlines()
        run_data.fasta_files = [GenomeData.from_path(entry, "nt-fasta-file") for entry in entries_list]
        for gd in run_data.fasta_files:
            gd.prodigal_used = True

    if args.amino_acid_files:
        with open(args.amino_acid_files, "r") as f:
            entries_list = f.read().splitlines()
        run_data.amino_acid_files = [GenomeData.from_path(entry, "aa-fasta-file") for entry in entries_list]

    run_data.update_all_input_genomes()
    run_data.run_files_dir = args.run_files_dir
    run_data.run_files_dir_rel = args.run_files_dir_rel
    run_data.run_data_path = run_data.run_files_dir + "/run-data.json"

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
    run_data_dict = asdict(run_data)
    if isinstance(run_data.start_time, datetime):
        run_data_dict['start_time'] = run_data.start_time.isoformat()
    with open(run_data.run_data_path, "w") as f:
        json.dump(run_data_dict, f, indent=2)


def read_run_data(path):
    try:
        with open(path, "r") as f:
            run_data_dict = json.load(f)

        if "ncbi_accs" in run_data_dict:
            run_data_dict["ncbi_accs"] = [GenomeData(**gd) if isinstance(gd, dict) else gd for gd in run_data_dict["ncbi_accs"]]
        if "genbank_files" in run_data_dict:
            run_data_dict["genbank_files"] = [GenomeData(**gd) if isinstance(gd, dict) else gd for gd in run_data_dict["genbank_files"]]
        if "fasta_files" in run_data_dict:
            run_data_dict["fasta_files"] = [GenomeData(**gd) if isinstance(gd, dict) else gd for gd in run_data_dict["fasta_files"]]
        if "amino_acid_files" in run_data_dict:
            run_data_dict["amino_acid_files"] = [GenomeData(**gd) if isinstance(gd, dict) else gd for gd in run_data_dict["amino_acid_files"]]
        if "all_input_genomes" in run_data_dict:
            run_data_dict["all_input_genomes"] = [GenomeData(**gd) if isinstance(gd, dict) else gd for gd in run_data_dict["all_input_genomes"]]

        if "tools_used" in run_data_dict and run_data_dict["tools_used"] is not None:
            run_data_dict["tools_used"] = ToolsUsed(**run_data_dict["tools_used"])
        else:
            run_data_dict["tools_used"] = ToolsUsed()

        if "start_time" in run_data_dict:
            run_data_dict["start_time"] = datetime.fromisoformat(run_data_dict["start_time"])

        if "SCG_targets" in run_data_dict:
            run_data_dict["SCG_targets"] = [SCGset(**scg) if isinstance(scg, dict) else scg for scg in run_data_dict["SCG_targets"]]
        run_data = RunData(**run_data_dict)
        run_data.update_all_input_genomes()

        return run_data

    except FileNotFoundError:
        pass


def get_snakefile_path(basename):
    return resource_filename("gtotree", f"smk/{basename}")


def touch(path):
    with open(path, 'a'):
        os.utime(path, None)


def run_snakemake(snakefile, tqdm_jobs, args, run_data, description, print_lines=False):
    print("")
    num_finished = 0 # counting this so it doesn't flip the progress bar when it counts the final "all" rule
    log = f"{run_data.logs_dir}/{description.replace(' ', '-').lower()}.log"
    log_dir_rel = f"{run_data.logs_dir_rel}/{description.replace(' ', '-').lower()}.log"

    cmd = [
        "snakemake",
        "--snakefile", snakefile,
        "--cores", f"{args.num_jobs}",
        "--default-resources", f"tmpdir='{run_data.tmp_dir}'",
        "--config", f"run_data_path={run_data.run_data_path}",
        "--directory", run_data.run_files_dir
    ]

    bar_format = "      {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]"
    with open(log, "w") as log_file, tqdm(total = tqdm_jobs, bar_format = bar_format, ncols = 76) as pbar:
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
            report_snakemake_failure(description, log_dir_rel, run_data)
            report_early_exit(run_data)


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


def check_file_exists_and_not_empty(path):
    try:
        if os.path.getsize(path) > 0:
            return True
        else:
            remove_file_if_exists(path)
    except FileNotFoundError:
        pass
    return False


def concat_files(file_list, output_file):
    with open(output_file, 'w') as outfile:
        for fname in file_list:
            with open(fname, 'r') as infile:
                shutil.copyfileobj(infile, outfile)
