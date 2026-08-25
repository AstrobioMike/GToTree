import os
import pytest  # type: ignore
from gtotree.utils.misc.general import GenomeData, RunData, apply_suffixes_to_mapping_dict
from gtotree.utils.misc.preflight_checks import (build_mapping_key_lookup,
                                            normalize_mapping_dict_keys,
                                            check_all_mapping_file_entries_are_in_input_genomes,
                                            make_mapping_dict)


def _run_data():
    rd = RunData()
    rd.mapping_file_path = "map.tsv"
    rd.amino_acid_files = [
        GenomeData.from_path("genomes/mock-1.faa", "amino-acid-fasta"),
        GenomeData.from_path("genomes/mock-2.faa", "amino-acid-fasta"),
    ]
    rd.ncbi_accs = [GenomeData.from_acc("GCF_000000001.1")]
    rd.update_all_input_genomes()
    return rd


class TestMappingKeyLookup:

    def test_accepts_the_path_as_provided(self):
        lookup = build_mapping_key_lookup(_run_data())
        assert lookup["genomes/mock-1.faa"] == "mock-1"

    def test_accepts_the_genome_id(self):
        lookup = build_mapping_key_lookup(_run_data())
        assert lookup["mock-1"] == "mock-1"

    def test_an_accession_is_its_own_id(self):
        lookup = build_mapping_key_lookup(_run_data())
        assert lookup["GCF_000000001.1"] == "GCF_000000001.1"

    def test_covers_every_input_genome(self):
        rd = _run_data()
        lookup = build_mapping_key_lookup(rd)
        assert set(lookup.values()) == {gd.id for gd in rd.all_input_genomes}


class TestNormalizeMappingDictKeys:

    def test_rekeys_paths_to_genome_ids(self):
        """
        The bug this fixes: alignment headers are genome ids, so a dict keyed by
        the user's paths never matched and -m silently did nothing for -g/-f/-A.
        """
        rd = _run_data()
        raw = {"genomes/mock-1.faa": "Alpha_one", "genomes/mock-2.faa": "Beta_two"}

        assert normalize_mapping_dict_keys(raw, rd) == {
            "mock-1": "Alpha_one",
            "mock-2": "Beta_two",
        }

    def test_accessions_are_unchanged(self):
        """Why this went unnoticed: for -a the two key forms are identical."""
        rd = _run_data()
        raw = {"GCF_000000001.1": "Some_label"}
        assert normalize_mapping_dict_keys(raw, rd) == raw

    def test_ids_written_directly_are_left_alone(self):
        rd = _run_data()
        raw = {"mock-1": "Alpha_one"}
        assert normalize_mapping_dict_keys(raw, rd) == {"mock-1": "Alpha_one"}

    def test_a_mixed_mapping_file_lands_in_one_key_space(self):
        rd = _run_data()
        raw = {"genomes/mock-1.faa": "Alpha_one",
               "mock-2": "Beta_two",
               "GCF_000000001.1": "Gamma_three"}

        normalized = normalize_mapping_dict_keys(raw, rd)

        assert set(normalized) <= {gd.id for gd in rd.all_input_genomes}
        assert normalized["mock-1"] == "Alpha_one"

    def test_key_order_is_preserved(self):
        rd = _run_data()
        raw = {"genomes/mock-2.faa": "Beta_two", "genomes/mock-1.faa": "Alpha_one"}
        assert list(normalize_mapping_dict_keys(raw, rd)) == ["mock-2", "mock-1"]

    def test_values_are_untouched(self):
        rd = _run_data()
        raw = {"genomes/mock-1.faa": "Alpha_one"}
        assert list(normalize_mapping_dict_keys(raw, rd).values()) == ["Alpha_one"]


class TestValidationUsesTheSameAcceptedKeys:
    """
    Validator and normalizer read the same lookup, so anything accepted is
    guaranteed resolvable -- they can't drift into accepting a form that then
    fails to rekey.
    """

    @pytest.mark.parametrize("key", ["genomes/mock-1.faa", "mock-1",
                                     "GCF_000000001.1"])
    def test_accepted_keys_pass_validation(self, key):
        rd = _run_data()
        check_all_mapping_file_entries_are_in_input_genomes({key: "label"}, rd)

    def test_an_unknown_entry_is_still_rejected(self):
        rd = _run_data()
        with pytest.raises(SystemExit):
            check_all_mapping_file_entries_are_in_input_genomes(
                {"not-a-genome.faa": "label"}, rd)

    @pytest.mark.parametrize("key", ["genomes/mock-1.faa", "mock-1",
                                     "GCF_000000001.1"])
    def test_everything_validation_accepts_can_be_rekeyed(self, key):
        rd = _run_data()
        normalized = normalize_mapping_dict_keys({key: "label"}, rd)
        assert set(normalized) <= {gd.id for gd in rd.all_input_genomes}


def _write(tmp_path, text):
    p = tmp_path / "mapping.tsv"
    p.write_text(text)
    return str(p)


def test_suffix_only_row_is_deferred_not_baked(tmp_path):
    # column 2 empty, column 3 present -> suffix-only, must not land in mapping_dict
    path = _write(tmp_path, "mike\tlee\t\nGCA_049325295.1\t\tsuffix\n")
    mapping_dict, suffix_dict = make_mapping_dict(path)

    assert mapping_dict == {"mike": "lee"}
    assert suffix_dict == {"GCA_049325295.1": "suffix"}
    # the whole point: the suffixed genome is absent from mapping_dict so the
    # taxonomy handler is free to build a base label for it
    assert "GCA_049325295.1" not in mapping_dict


def test_full_replacement_with_suffix_is_still_baked(tmp_path):
    # column 2 present alongside column 3 -> full label, suffix baked in now
    path = _write(tmp_path, "GCA_1\tGenus species\ttag\n")
    mapping_dict, suffix_dict = make_mapping_dict(path)

    assert mapping_dict == {"GCA_1": "Genus_species_tag"}
    assert suffix_dict == {}


class _RD:
    def __init__(self, mapping_dict, suffix_dict):
        self.mapping_dict = dict(mapping_dict)
        self.suffix_dict = dict(suffix_dict)


def test_suffix_appends_onto_taxonomy_base():
    # base label was built by the taxonomy handler; suffix goes on the end
    rd = _RD(
        {"GCA_049325295.1": "GCA_049325295.1_Bacteria_Escherichia_coli"},
        {"GCA_049325295.1": "suffix"},
    )
    apply_suffixes_to_mapping_dict(rd)
    assert rd.mapping_dict["GCA_049325295.1"] == \
        "GCA_049325295.1_Bacteria_Escherichia_coli_suffix"


def test_suffix_appends_onto_bare_id_when_no_taxonomy():
    # no taxonomy label was built; suffix goes onto the genome id itself
    rd = _RD({}, {"GCA_049325295.1": "suffix"})
    apply_suffixes_to_mapping_dict(rd)
    assert rd.mapping_dict["GCA_049325295.1"] == "GCA_049325295.1_suffix"


def test_apply_suffixes_is_a_noop_without_suffix_entries():
    rd = _RD({"GCA_1": "Whatever"}, {})
    apply_suffixes_to_mapping_dict(rd)
    assert rd.mapping_dict == {"GCA_1": "Whatever"}
