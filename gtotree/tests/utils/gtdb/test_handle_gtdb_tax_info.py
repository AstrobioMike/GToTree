import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pytest # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.gtdb.get_gtdb_data import PARQUET_FILENAME
import gtotree.utils.gtdb.handle_gtdb_tax_info as mod


def _write_mock_gtdb(path, records):
    """records: list of (genbank_acc, domain..species tuple). Builds a schema-faithful
    slim GTDB parquet with just the columns the annotator reads."""
    cols = {"ncbi_genbank_assembly_accession": []}
    for r in RANKS:
        cols[r] = []
    for gb_acc, ranks in records:
        cols["ncbi_genbank_assembly_accession"].append(gb_acc)
        for i, r in enumerate(RANKS):
            cols[r].append(ranks[i])
    pq.write_table(pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}),
                   str(path))


class _RunData:
    def __init__(self, accs):
        self._accs = accs
        self.mapping_dict = {}
        self.tmp_dir = "/tmp"

    def remaining_ncbi_accs(self):
        return self._accs


class _Args:
    lineage = "Domain,Phylum,Genus,Species"


@pytest.fixture
def gtdb_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    _write_mock_gtdb(
        tmp_path / PARQUET_FILENAME,
        [
            ("GCA_000005845.2",
             ("Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales",
              "Enterobacteriaceae", "Escherichia", "Escherichia coli")),
            ("GCA_000009065.1",
             ("Bacteria", "Bacillota", "Bacilli", "Bacillales",
              "Bacillaceae", "Bacillus", "Bacillus subtilis")),
        ],
    )
    return str(tmp_path)


def test_gca_input_is_annotated(gtdb_dir):
    rd = _RunData(["GCA_000005845.2"])
    rd = mod.update_mapping_dict_with_gtdb_tax_info(_Args(), rd)
    # Domain,Phylum,Genus,Species with genus-stripped species
    assert rd.mapping_dict["GCA_000005845.2"] == \
        "GCA_000005845.2_Bacteria_Pseudomonadota_Escherichia_coli"


def test_gcf_input_matches_gtdb_gca_record(gtdb_dir):
    # GTDB stores GCA_000005845.2; an input GCF_ with the same numeric core must match
    rd = _RunData(["GCF_000005845.2"])
    rd = mod.update_mapping_dict_with_gtdb_tax_info(_Args(), rd)
    assert rd.mapping_dict["GCF_000005845.2"] == \
        "GCF_000005845.2_Bacteria_Pseudomonadota_Escherichia_coli"


def test_species_genus_prefix_is_stripped(gtdb_dir):
    rows = mod.subset_gtdb_info(_RunData(["GCA_000009065.1"]))
    assert rows["GCA_000009065.1"]["species"] == "subtilis"
    assert rows["GCA_000009065.1"]["genus"] == "Bacillus"


def test_unmatched_accession_is_absent(gtdb_dir):
    rd = _RunData(["GCA_999999999.9"])
    rd = mod.update_mapping_dict_with_gtdb_tax_info(_Args(), rd)
    assert rd.mapping_dict == {}


def test_no_input_accessions_returns_empty(gtdb_dir):
    assert mod.subset_gtdb_info(_RunData([])) == {}


def test_existing_mapping_is_not_overwritten(gtdb_dir):
    rd = _RunData(["GCA_000005845.2"])
    rd.mapping_dict["GCA_000005845.2"] = "already_here"
    rd = mod.update_mapping_dict_with_gtdb_tax_info(_Args(), rd)
    assert rd.mapping_dict["GCA_000005845.2"] == "already_here"


def test_strain_rank_is_skipped_in_label(gtdb_dir):
    class ArgsWithStrain:
        lineage = "Domain,Strain,Species"
    rd = _RunData(["GCA_000005845.2"])
    rd = mod.update_mapping_dict_with_gtdb_tax_info(ArgsWithStrain(), rd)
    # 'strain' has no column in GTDB and is skipped -> Domain + Species only
    assert rd.mapping_dict["GCA_000005845.2"] == "GCA_000005845.2_Bacteria_coli"
