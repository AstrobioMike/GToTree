import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pyarrow.compute as pc # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS, accession_core
from gtotree.utils.gtdb.get_gtdb_data import gtdb_data_table_path

# the 7 rank columns plus the join key, read from the Parquet asset
_WANTED_COLUMNS = ["ncbi_genbank_assembly_accession"] + list(RANKS)

# extracts the numeric core (the run of digits) from a full assembly accession like
    # "GCA_000739065.1" -> "000739065"
_ACCESSION_CORE_REGEX = r'^[A-Za-z]+_([0-9]+)\.[0-9]+$'


def update_mapping_dict_with_gtdb_tax_info(args, run_data):

    tax_rows = subset_gtdb_info(run_data)

    run_data.tax_info_dict = tax_rows

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
    """
    lineages = gtdb_lineages_for_accessions(run_data.remaining_ncbi_accs())

    for ranks in lineages.values():
        # removing 'genus' from 'species' otherwise if the user wants genus and species,
        # it would list the genus twice
        genus = ranks["genus"]
        species = ranks["species"]
        if isinstance(species, str):
            ranks["species"] = species.removeprefix(f"{genus} ")

    return lineages


def gtdb_lineages_for_accessions(accessions, gtdb_path=None):
    """
    Map assembly accessions to their GTDB lineages, verbatim from the Parquet asset.

    Accessions with no GTDB row are simply absent from the result

    Returns a dict: input_acc -> {rank: value for the 7 ranks}.

    Shared by the two callers that need this join:
    header re-labelling (`-D`) and the `-H` auto-pick, which links an NCBI-sourced `-w`
    selection back to the GTDB taxonomy the pre-built SCG-sets are built from
    """
    # map each wanted accession to its numeric core for the join
    core_to_input = {}
    for acc in accessions:
        core = accession_core(acc)
        if core:
            core_to_input[core] = acc

    if not core_to_input:
        return {}

    if gtdb_path is None:
        gtdb_path = gtdb_data_table_path()

    wanted_cores = pa.array(list(core_to_input.keys()))

    result = {}
    parquet_file = pq.ParquetFile(gtdb_path)

    for group_index in range(parquet_file.metadata.num_row_groups):
        group = parquet_file.read_row_group(group_index, columns=_WANTED_COLUMNS)

        cores = pc.replace_substring_regex(
            group.column("ncbi_genbank_assembly_accession"),
            _ACCESSION_CORE_REGEX, r'\1')
        matched = group.filter(pc.is_in(cores, value_set=wanted_cores))

        if matched.num_rows == 0:
            continue

        gb_matched = matched.column("ncbi_genbank_assembly_accession").to_pylist()
        rank_cols = {rank: matched.column(rank).to_pylist() for rank in RANKS}

        for i, gb_acc in enumerate(gb_matched):
            input_acc = core_to_input.get(accession_core(gb_acc))
            if input_acc is None:
                continue
            result[input_acc] = {rank: rank_cols[rank][i] for rank in RANKS}

    return result
