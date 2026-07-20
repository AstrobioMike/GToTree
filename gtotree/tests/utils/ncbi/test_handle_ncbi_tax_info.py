import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pytest # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.ncbi.get_ncbi_assembly_data import PARQUET_FILENAME
import gtotree.utils.ncbi.handle_ncbi_tax_info as mod


def _write_mock_ncbi(path, records):
    """records: list of (assembly_accession, lineage-7-tuple, infraspecific_name)."""
    cols = {"assembly_accession": [], "infraspecific_name": []}
    for r in RANKS:
        cols[r] = []
    for acc, lineage, infra in records:
        cols["assembly_accession"].append(acc)
        cols["infraspecific_name"].append(infra)
        for i, r in enumerate(RANKS):
            cols[r].append(lineage[i])
    pq.write_table(pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}),
                   str(path))


class _RunData:
    def __init__(self, accs):
        self._accs = accs
        self.mapping_dict = {}

    def remaining_ncbi_accs(self):
        return self._accs


class _Args:
    lineage = "Domain,Phylum,Genus,Species"


class _ArgsStrain:
    lineage = "Domain,Genus,Species,Strain"


@pytest.fixture
def ncbi_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
    _write_mock_ncbi(
        tmp_path / PARQUET_FILENAME,
        [
            ("GCF_000005845.2",
             ("Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales",
              "Enterobacteriaceae", "Escherichia", "Escherichia coli"),
             "strain=K-12"),
            ("GCA_000009065.1",
             ("Bacteria", "Bacillota", "Bacilli", "Bacillales",
              "Bacillaceae", "Bacillus", "Bacillus subtilis"),
             "na"),
        ],
    )
    return str(tmp_path)


def test_basic_annotation(ncbi_dir):
    rd = _RunData(["GCF_000005845.2"])
    rd = mod.update_mapping_dict_with_ncbi_tax_info(_Args(), rd)
    assert rd.mapping_dict["GCF_000005845.2"] == \
        "GCF_000005845.2_Bacteria_Pseudomonadota_Escherichia_coli"


def test_gca_input_matches_gcf_row_by_core(ncbi_dir):
    # input GCA matches the GCF row via accession core
    rd = _RunData(["GCA_000005845.2"])
    rd = mod.update_mapping_dict_with_ncbi_tax_info(_Args(), rd)
    assert rd.mapping_dict["GCA_000005845.2"] == \
        "GCA_000005845.2_Bacteria_Pseudomonadota_Escherichia_coli"


def test_species_genus_stripped(ncbi_dir):
    rows = mod.subset_ncbi_info(_RunData(["GCF_000005845.2"]))
    assert rows["GCF_000005845.2"]["species"] == "coli"
    assert rows["GCF_000005845.2"]["genus"] == "Escherichia"


def test_strain_from_infraspecific_name(ncbi_dir):
    rd = _RunData(["GCF_000005845.2"])
    rd = mod.update_mapping_dict_with_ncbi_tax_info(_ArgsStrain(), rd)
    # Domain,Genus,Species,Strain -> strain K-12 appended
    assert rd.mapping_dict["GCF_000005845.2"] == \
        "GCF_000005845.2_Bacteria_Escherichia_coli_K-12"


def test_strain_na_when_infraspecific_is_na(ncbi_dir):
    rows = mod.subset_ncbi_info(_RunData(["GCA_000009065.1"]))
    assert rows["GCA_000009065.1"]["strain"] == "NA"


def test_strain_spaces_collapsed_to_hyphens(tmp_path, monkeypatch):
    monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
    _write_mock_ncbi(
        tmp_path / PARQUET_FILENAME,
        [("GCF_111111111.1",
          ("Bacteria", "P", "C", "O", "F", "Genusy", "Genusy speciesy"),
          "strain=SCGC AAA011-O16")],
    )
    rows = mod.subset_ncbi_info(_RunData(["GCF_111111111.1"]))
    assert rows["GCF_111111111.1"]["strain"] == "SCGC-AAA011-O16"


def test_unmatched_returns_empty(ncbi_dir):
    rd = _RunData(["GCF_999999999.9"])
    rd = mod.update_mapping_dict_with_ncbi_tax_info(_Args(), rd)
    assert rd.mapping_dict == {}


def test_no_input_returns_empty(ncbi_dir):
    assert mod.subset_ncbi_info(_RunData([])) == {}


def test_existing_mapping_not_overwritten(ncbi_dir):
    rd = _RunData(["GCF_000005845.2"])
    rd.mapping_dict["GCF_000005845.2"] = "already_here"
    rd = mod.update_mapping_dict_with_ncbi_tax_info(_Args(), rd)
    assert rd.mapping_dict["GCF_000005845.2"] == "already_here"
