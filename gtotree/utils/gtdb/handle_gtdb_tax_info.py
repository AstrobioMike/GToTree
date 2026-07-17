import os
import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pyarrow.compute as pc # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS, accession_core
from gtotree.utils.gtdb.get_gtdb_data import gtdb_data_table_path

# the 7 rank columns plus the join key, read from the Parquet asset
_WANTED_COLUMNS = ["ncbi_genbank_assembly_accession"] + list(RANKS)

# extracts the numeric core (the run of digits) from a full assembly accession like
# "GCA_000739065.1" -> "000739065"; used Arrow-side to build the join key column
_ACCESSION_CORE_REGEX = r'^[A-Za-z]+_([0-9]+)\.[0-9]+$'


def update_mapping_dict_with_gtdb_tax_info(args, run_data):

    tax_rows = subset_gtdb_info(run_data)

    for input_acc, ranks in tax_rows.items():
        if input_acc in run_data.mapping_dict:
            continue

        new_label = f"{input_acc}_"
        for rank in args.lineage.split(","):
            if rank.lower() != "strain":
                new_label += ranks[rank.lower()] + "_"

        new_label = new_label[:-1]           # removing last underscore
        new_label = new_label.replace(" ", "_")
        run_data.mapping_dict[input_acc] = new_label

    return run_data


def subset_gtdb_info(run_data):
    """
    Pull GTDB lineage info for the remaining NCBI input accessions from the hosted
    Parquet asset.

    Returns a dict: input_acc -> {rank: value for the 7 ranks}, with the species value
    de-duplicated of its leading genus

    The join is done on the version-stripped numeric accession core
    """
    input_accs = run_data.remaining_ncbi_accs()

    # map each wanted accession to its numeric core for the join
    core_to_input = {}
    for acc in input_accs:
        core = accession_core(acc)
        if core:
            core_to_input[core] = acc

    if not core_to_input:
        return {}

    gtdb_path = gtdb_data_table_path(os.environ['GTDB_DIR'])

    gb_col = pq.read_table(gtdb_path, columns=["ncbi_genbank_assembly_accession"]).column(0)
    cores = pc.replace_substring_regex(gb_col, _ACCESSION_CORE_REGEX, r'\1')
    mask = pc.is_in(cores, value_set=pa.array(list(core_to_input.keys())))
    keep_idx = pc.indices_nonzero(mask)

    subset = pq.read_table(gtdb_path, columns=_WANTED_COLUMNS).take(keep_idx)

    gb_matched = subset.column("ncbi_genbank_assembly_accession").to_pylist()
    rank_cols = {rank: subset.column(rank).to_pylist() for rank in RANKS}

    result = {}
    for i, gb_acc in enumerate(gb_matched):
        input_acc = core_to_input.get(accession_core(gb_acc))
        if input_acc is None:
            continue
        ranks = {rank: rank_cols[rank][i] for rank in RANKS}
        # removing 'genus' from 'species' otherwise if the user wants genus and species,
        # it would list the genus twice
        genus = ranks["genus"]
        species = ranks["species"]
        if isinstance(species, str):
            ranks["species"] = species.removeprefix(f"{genus} ")
        result[input_acc] = ranks

    return result
