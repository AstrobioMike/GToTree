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


def test_bacteria_is_leftmost(small_table):
    _, domains = G.group_sets(G.parse_sets(small_table))
    assert [d[0] for d in domains] == ["Bacteria", "Archaea"]


def test_domain_order_holds_when_the_domains_are_the_same_size():
    """Mid-build the two can tie, and an alphabetical tie-break would put Bacteria on
    the right with nothing to its left."""
    rows = [entry("Bdellovibrionota"), entry("Asgardarchaeota", domain="Archaea")]
    _, domains = G.group_sets(G.parse_sets(make_df(rows)))
    assert [d[0] for d in domains] == ["Bacteria", "Archaea"]


def test_an_unrecognized_domain_sorts_after_the_known_ones_by_size():
    rows = [entry("Bdellovibrionota"), entry("Asgardarchaeota", domain="Archaea"),
            entry("Somethingelse", domain="Eukarya")]
    _, domains = G.group_sets(G.parse_sets(make_df(rows)))
    assert [d[0] for d in domains] == ["Bacteria", "Archaea", "Eukarya"]


def test_the_domain_level_set_becomes_the_heading(small_table):
    _, domains = G.group_sets(G.parse_sets(small_table))
    _, domain_set, blocks = domains[0]
    assert domain_set["name"] == "Bacteria"
    assert "Bacteria" not in [parent["name"] for parent, _ in blocks]


def test_a_class_nests_under_its_parent_phylum(small_table):
    _, domains = G.group_sets(G.parse_sets(small_table))
    blocks = dict((parent["name"], [k["name"] for k in kids])
                  for parent, kids in domains[0][2])
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


def test_the_domains_sit_side_by_side():
    lines = G.layout_lines([], domains_for(36))
    heading = next(line for line in lines if "Bacteria" in line)
    assert "Archaea" in heading


def test_every_set_appears_exactly_once():
    text = "\n".join(G.layout_lines([], domains_for(36)))
    for i in range(36):
        assert text.count(f"Bacterialphylum{i:02d}") == 1


def test_each_domain_gets_one_entry_per_line():
    """No reflowing into extra columns -- one column per domain, always."""
    lines = G.layout_lines([], domains_for(36))
    assert all(line.count("Bacterialphylum") <= 1 for line in lines)


def test_the_listing_is_as_tall_as_the_biggest_domain():
    domains = domains_for(36, n_archaeal_phyla=7)
    lines = [line for line in G.layout_lines([], domains) if line.strip()]
    assert len(lines) == 37  # the Bacteria heading plus its 36 phyla


def test_a_domain_that_runs_out_leaves_the_short_side_on_the_right():
    """Bacteria stays at the margin once Archaea has run out of entries, rather than
    the surviving column being stranded off to the right."""
    lines = G.layout_lines([], domains_for(36, n_archaeal_phyla=2))
    tail = [line for line in lines if "Bacterialphylum" in line][-1]
    assert "Archaealphylum" not in tail
    assert tail.index("Bacterialphylum") == G.LEFT_MARGIN + 2


def test_columns_stay_aligned_when_a_domain_runs_out_of_entries():
    lines = G.layout_lines([], domains_for(36, n_archaeal_phyla=2))
    body = [line for line in lines if "Archaealphylum" in line]
    assert len({line.index("Archaealphylum") for line in body}) == 1


def test_cross_domain_sets_are_rendered_above_the_domains(small_table):
    cross, domains = G.group_sets(G.parse_sets(small_table))
    text = "\n".join(G.layout_lines(cross, domains))
    assert text.index("Universal-Hug-et-al") < text.index("Halobacteriota")


def test_a_domain_with_no_domain_level_set_still_gets_a_heading():
    rows = [entry("Bacteroidota"), entry("Zixibacteria")]
    _, domains = G.group_sets(G.parse_sets(make_df(rows)))
    assert domains[0][1] is None
    assert any("Bacteria" in line for line in G.layout_lines([], domains))


def test_layout_is_empty_when_there_is_nothing_to_lay_out():
    assert G.layout_lines([], []) == []


def test_no_line_has_trailing_whitespace(small_table):
    cross, domains = G.group_sets(G.parse_sets(small_table))
    assert all(line == line.rstrip() for line in G.layout_lines(cross, domains))


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


################################################################################
# alignment
################################################################################

def count_columns(lines):
    """Where each `(n)` sits on each line."""
    positions = set()
    for line in lines:
        start = 0
        while True:
            start = line.find("(", start)
            if start < 0:
                break
            end = line.find(")", start)
            positions.add(end)
            start = end + 1
    return positions


def test_counts_line_up_across_the_whole_listing(small_table):
    """A count in the top block that doesn't line up with a count in a domain strip
    below it reads as a mistake, so everything is measured once, together."""
    cross, domains = G.group_sets(G.parse_sets(small_table))
    lines = G.layout_lines(cross, domains)
    positions = sorted(count_columns(lines))
    # one position per column of the widest strip, and the top block shares the first
    assert positions[0] == min(positions)
    spacings = {b - a for a, b in zip(positions, positions[1:])}
    assert len(spacings) <= 1


def test_the_top_block_shares_the_left_margin_with_the_domains(small_table):
    cross, domains = G.group_sets(G.parse_sets(small_table))
    lines = [line for line in G.layout_lines(cross, domains) if line.strip()]
    indents = {len(line) - len(line.lstrip()) for line in lines}
    assert min(indents) == G.LEFT_MARGIN


def test_domain_headings_are_not_singled_out_with_colour(small_table):
    """They're ordinary selectable sets once they're built, so nothing should imply
    otherwise."""
    cross, domains = G.group_sets(G.parse_sets(small_table))
    assert not any("\033[" in line for line in G.layout_lines(cross, domains))
