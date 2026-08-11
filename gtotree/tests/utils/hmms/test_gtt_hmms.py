import pandas as pd  # type: ignore
import pytest

from gtotree.utils.hmms import gtt_hmms as G


V1_COLUMNS = ["file", "source", "target_taxa", "num_genes", "link", "md5"]


def make_df(rows, columns=None):
    """A stand-in info table. `columns` forces a schema, e.g. the v1 one."""
    frame = pd.DataFrame(rows)
    if columns is not None:
        for column in columns:
            if column not in frame.columns:
                frame[column] = "x"
        frame = frame[columns]
    return frame


def entry(name, num_genes=10, rank="phylum", domain="Bacteria", parent="NA",
          gtdb_release="r232"):
    return {"file": f"{name}.hmm", "num_genes": num_genes, "rank": rank,
            "domain": domain, "parent": parent, "gtdb_release": gtdb_release}


@pytest.fixture
def small_table():
    return make_df([
        entry("Universal-Hug-et-al", 16, "universal", "NA"),
        entry("Bacteria-and-Archaea", 25, "multi-domain", "NA"),
        entry("Bacteria", 74, "domain", "Bacteria"),
        entry("Archaea", 76, "domain", "Archaea"),
        entry("Pseudomonadota", 119),
        entry("Gammaproteobacteria", 172, "class", "Bacteria", "Pseudomonadota"),
        entry("Bacteroidota", 90),
        entry("Halobacteriota", 60, "phylum", "Archaea"),
    ])


################################################################################
# reading the table
################################################################################

def test_has_layout_info_true_for_a_v2_table(small_table):
    assert G.has_layout_info(small_table)


def test_has_layout_info_false_for_the_v1_schema(small_table):
    assert not G.has_layout_info(make_df(small_table.to_dict("records"), V1_COLUMNS))


def test_has_layout_info_false_when_the_columns_are_there_but_empty():
    """A table with the columns added but never populated should list flat, not as
    one unlabelled heap."""
    rows = [entry("Bacteroidota", rank="NA", domain="NA")]
    assert not G.has_layout_info(make_df(rows))


def test_parse_sets_strips_the_suffix_and_normalizes_absent_values():
    parsed = G.parse_sets(make_df([entry("Bacteroidota", parent="NA")]))
    assert parsed[0]["name"] == "Bacteroidota"
    assert parsed[0]["parent"] == ""


def test_parse_sets_treats_a_missing_gene_count_as_unknown():
    parsed = G.parse_sets(make_df([entry("Bacteroidota", num_genes="NA")]))
    assert parsed[0]["num_genes"] == "?"


################################################################################
# grouping
################################################################################

def test_cross_domain_sets_lead_broadest_first(small_table):
    cross, _ = G.group_sets(G.parse_sets(small_table))
    assert [s["name"] for s in cross] == ["Universal-Hug-et-al", "Bacteria-and-Archaea"]


def test_domains_are_ordered_smallest_first(small_table):
    _, domains = G.group_sets(G.parse_sets(small_table))
    assert [d[0] for d in domains] == ["Archaea", "Bacteria"]


def test_the_domain_level_set_becomes_the_heading(small_table):
    _, domains = G.group_sets(G.parse_sets(small_table))
    _, domain_set, blocks = domains[1]
    assert domain_set["name"] == "Bacteria"
    assert "Bacteria" not in [parent["name"] for parent, _ in blocks]


def test_a_class_nests_under_its_parent_phylum(small_table):
    _, domains = G.group_sets(G.parse_sets(small_table))
    blocks = dict((parent["name"], [k["name"] for k in kids])
                  for parent, kids in domains[1][2])
    assert blocks["Pseudomonadota"] == ["Gammaproteobacteria"]
    assert blocks["Bacteroidota"] == []


def test_a_class_whose_parent_was_not_built_still_gets_listed():
    """Nothing may go missing just because its parent didn't make the cut."""
    rows = [entry("Bacteria", 74, "domain", "Bacteria"),
            entry("Gammaproteobacteria", 172, "class", "Bacteria", "Pseudomonadota")]
    _, domains = G.group_sets(G.parse_sets(make_df(rows)))
    parents = [parent["name"] for parent, _ in domains[0][2]]
    assert parents == ["Gammaproteobacteria"]


def test_phyla_within_a_domain_are_alphabetical():
    rows = [entry(name) for name in ("Zixibacteria", "Acidobacteriota", "Myxococcota")]
    _, domains = G.group_sets(G.parse_sets(make_df(rows)))
    assert [parent["name"] for parent, _ in domains[0][2]] == [
        "Acidobacteriota", "Myxococcota", "Zixibacteria"]


def test_a_domain_with_no_domain_level_set_still_renders(small_table):
    rows = [entry("Bacteroidota"), entry("Zixibacteria")]
    _, domains = G.group_sets(G.parse_sets(make_df(rows)))
    assert domains[0][1] is None
    assert any("Bacteria" in line for line in G.layout_lines([], domains, 100))


################################################################################
# column packing
################################################################################

def test_pack_blocks_never_splits_a_block_across_columns():
    blocks = [["a1", "a2", "a3"], ["b1"], ["c1", "c2"], ["d1"], ["e1"]]
    columns = G.pack_blocks(blocks, 2)
    joined = ["".join(column) for column in columns]
    for block in blocks:
        assert sum("".join(block) in column for column in joined) == 1


def test_pack_blocks_preserves_order_reading_down_then_across():
    blocks = [[str(i)] for i in range(6)]
    assert G.pack_blocks(blocks, 2) == [["0", "1", "2"], ["3", "4", "5"]]


def test_pack_blocks_always_returns_the_requested_number_of_columns():
    assert len(G.pack_blocks([["a"]], 3)) == 3


################################################################################
# layout
################################################################################

def domains_for(n_bacterial_phyla, n_archaeal_phyla=7):
    rows = [entry("Bacteria", 74, "domain", "Bacteria"),
            entry("Archaea", 76, "domain", "Archaea")]
    rows += [entry(f"Bacterialphylum{i:02d}") for i in range(n_bacterial_phyla)]
    rows += [entry(f"Archaealphylum{i:02d}", domain="Archaea")
             for i in range(n_archaeal_phyla)]
    _, domains = G.group_sets(G.parse_sets(make_df(rows)))
    return domains


def test_a_wide_terminal_puts_the_domains_side_by_side():
    lines = G.layout_lines([], domains_for(36), 160)
    heading = next(line for line in lines if "Archaea" in line)
    assert "Bacteria" in heading


def test_every_set_appears_exactly_once_however_it_is_arranged():
    domains = domains_for(36)
    names = [f"Bacterialphylum{i:02d}" for i in range(36)]
    for width in (60, 80, 100, 120, 160, 200):
        text = "\n".join(G.layout_lines([], domains_for(36), width))
        for name in names:
            assert text.count(name) == 1, f"{name} at width {width}"


def test_the_chosen_arrangement_is_never_the_taller_one():
    """The point of the layout is to stay on screen, so whichever arrangement is
    shorter has to win -- including on a narrow terminal, where side by side leaves
    Bacteria stuck in a single tall column."""
    for width in (60, 80, 100, 120, 160):
        domains = domains_for(36)
        chosen = len(G.layout_lines([], domains, width))
        side = G._side_by_side([G.build_section(*d) for d in domains], width)
        stacked = G._stacked([G.build_section(*d) for d in domains], width)
        alternatives = [len(stacked)] + ([len(side)] if side is not None else [])
        assert chosen <= min(alternatives)


def test_a_narrow_terminal_stacks_the_domains():
    lines = G.layout_lines([], domains_for(36), 70)
    heading = next(line for line in lines if "Archaea" in line)
    assert "Bacteria" not in heading


def test_columns_stay_aligned_when_a_domain_runs_out_of_entries():
    lines = G.layout_lines([], domains_for(36, n_archaeal_phyla=2), 160)
    body = [line for line in lines if "Bacterialphylum" in line]
    starts = {line.index("Bacterialphylum") for line in body}
    assert len(starts) <= 2  # one per column of the Bacteria strip


def test_a_short_list_is_not_split_into_columns():
    """Four sets side by side is harder to read than four stacked, so a domain only
    earns a second column once it has enough entries to fill one."""
    lines = G.layout_lines([], domains_for(4, n_archaeal_phyla=3), 200)
    assert all(line.count("Bacterialphylum") <= 1 for line in lines)


def test_cross_domain_sets_are_rendered_above_the_domains(small_table):
    cross, domains = G.group_sets(G.parse_sets(small_table))
    lines = G.layout_lines(cross, domains, 160)
    text = "\n".join(lines)
    assert text.index("Universal-Hug-et-al") < text.index("Halobacteriota")


def test_layout_is_empty_when_there_is_nothing_to_lay_out():
    assert G.layout_lines([], [], 80) == []


################################################################################
# the flat fallback
################################################################################

def test_flat_listing_covers_every_set(capsys, small_table, tmp_path):
    G.report_available_scg_sets(str(tmp_path), small_table)
    out = capsys.readouterr().out
    for name in small_table["file"]:
        assert name.replace(".hmm", "") in out


def test_flat_listing_does_not_need_the_layout_columns(capsys, small_table, tmp_path):
    v1 = make_df(small_table.to_dict("records"), V1_COLUMNS)
    G.report_available_scg_sets(str(tmp_path), v1)
    assert "Pseudomonadota" in capsys.readouterr().out
