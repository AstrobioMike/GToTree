"""
Guards for `gtt itol` and the iToL writers it shares with post-search generation.

The writers are the kind of code whose failures are invisible from here: iToL
ignores a malformed row rather than reporting it, so a file that's wrong renders
as a tree with nothing on it. These tests pin the field layouts against the
published templates instead.
"""

import argparse
import os

import pytest  # type: ignore

from gtotree.cli.itol import build_parser, check_args
from gtotree.utils.misc.itol import (
    COLOR_MAP,
    SHAPE_MAP,
    ItolError,
    read_targets,
    resolve_color,
    write_binary_dataset,
    write_colorstrip,
    write_style_dataset,
    write_text_dataset,
)


@pytest.fixture
def targets_file(tmp_path):
    f = tmp_path / "targets.txt"
    f.write_text("genomeA\ngenomeB\n  genomeC  \n")
    return f


def _data_rows(path):
    lines = path.read_text().splitlines()
    return [l for l in lines[lines.index("DATA") + 1:] if l.strip()]


################################################################################
# read_targets
################################################################################

def test_read_targets_strips_blanks_and_duplicates(tmp_path):
    """
    A blank line becomes a row with an empty first field, which iToL drops along
    with whatever annotation was on it; a repeated ID becomes a duplicate row.
    """
    f = tmp_path / "messy.txt"
    f.write_text("genomeA\n\ngenomeB\ngenomeA\n  genomeC  \n\n")
    assert read_targets(str(f)) == ["genomeA", "genomeB", "genomeC"]


def test_read_targets_preserves_order(tmp_path):
    f = tmp_path / "ordered.txt"
    f.write_text("zeta\nalpha\nmu\n")
    assert read_targets(str(f)) == ["zeta", "alpha", "mu"]


def test_read_targets_missing_file_raises(tmp_path):
    with pytest.raises(ItolError, match="not found"):
        read_targets(str(tmp_path / "nope.txt"))


def test_read_targets_empty_file_raises(tmp_path):
    f = tmp_path / "empty.txt"
    f.write_text("\n  \n")
    with pytest.raises(ItolError, match="no genome IDs"):
        read_targets(str(f))


################################################################################
# resolve_color
################################################################################

@pytest.mark.parametrize("name", sorted(COLOR_MAP))
def test_resolve_color_accepts_every_named_color(name):
    assert resolve_color(name) == COLOR_MAP[name]


@pytest.mark.parametrize("hex_code", ["#0000ff", "#abc"])
def test_resolve_color_accepts_hex(hex_code):
    assert resolve_color(hex_code) == hex_code


@pytest.mark.parametrize("bad", ["chartreuse", "0000ff", "#00ff"])
def test_resolve_color_rejects_anything_else(bad):
    with pytest.raises(ItolError):
        resolve_color(bad)


################################################################################
# write_binary_dataset
################################################################################

@pytest.mark.parametrize("shape", sorted(SHAPE_MAP))
def test_binary_data_values_are_presence_flags_not_shape_codes(shape, tmp_path):
    """
    The DATASET_BINARY template allows only 1, 0 or -1 under DATA -- the symbol is
    chosen once, up in FIELD_SHAPES. Writing the shape code per row renders fine
    for `square` (code 1, which is also a valid presence flag) and silently draws
    nothing for every other shape, so this covers all of them.
    """
    out = tmp_path / "binary.txt"
    write_binary_dataset(str(out), ["gA", "gB"], shape=shape)

    text = out.read_text()
    assert text.startswith("DATASET_BINARY\nSEPARATOR TAB\n")
    assert f"FIELD_SHAPES\t{SHAPE_MAP[shape]}\n" in text
    assert _data_rows(out) == ["gA\t1", "gB\t1"]


def test_binary_field_labels_follow_the_dataset_label(tmp_path):
    """FIELD_LABELS is the text drawn above the symbol column."""
    out = tmp_path / "binary.txt"
    write_binary_dataset(str(out), ["gA"], label="my gene")
    text = out.read_text()
    assert "DATASET_LABEL\tmy gene\n" in text
    assert "FIELD_LABELS\tmy gene\n" in text


def test_binary_bad_shape_raises(tmp_path):
    with pytest.raises(ItolError, match="shapes are"):
        write_binary_dataset(str(tmp_path / "o.txt"), ["gA"], shape="hexagon")


################################################################################
# write_colorstrip
################################################################################

def test_colorstrip(tmp_path):
    out = tmp_path / "strip.txt"
    write_colorstrip(str(out), ["gA", "gB"], label="mygroup",
                     color=COLOR_MAP["green"], strip_width=40,
                     color_branches_too=True)

    text = out.read_text()
    assert text.startswith("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
    assert "COLOR_BRANCHES\t1\n" in text
    assert "STRIP_WIDTH\t40\n" in text
    assert _data_rows(out) == [f"{g}\t{COLOR_MAP['green']}\tmygroup"
                               for g in ("gA", "gB")]


def test_colorstrip_branches_off_by_default(tmp_path):
    out = tmp_path / "strip.txt"
    write_colorstrip(str(out), ["gA"])
    assert "COLOR_BRANCHES\t0\n" in out.read_text()


################################################################################
# write_style_dataset
################################################################################

def test_style_dataset_branches(tmp_path):
    out = tmp_path / "style.txt"
    write_style_dataset(str(out), ["gA", "gB"], label="my gene", color="#0000ff")

    text = out.read_text()
    assert text.startswith("DATASET_STYLE\nSEPARATOR TAB\n")
    assert "DATASET_LABEL\tmy gene\n" in text
    assert _data_rows(out) == [f"{g}\tbranch\tnode\t#0000ff\t3\tnormal"
                               for g in ("gA", "gB")]


def test_style_dataset_both_emits_branch_and_label_rows(tmp_path):
    out = tmp_path / "style.txt"
    write_style_dataset(str(out), ["gA", "gB"], what_to_color="both")
    rows = _data_rows(out)
    assert len(rows) == 4
    assert sum("\tbranch\t" in r for r in rows) == 2
    assert sum("\tlabel\t" in r for r in rows) == 2


################################################################################
# write_text_dataset
################################################################################

def test_text_dataset(tmp_path):
    out = tmp_path / "text.txt"
    write_text_dataset(str(out), ["gA", "gB"], "from site A", label="notes",
                       color="#0000ff")

    text = out.read_text()
    assert text.startswith("DATASET_TEXT\nSEPARATOR TAB\n")
    assert "DATASET_LABEL\tnotes\n" in text
    assert _data_rows(out) == [f"{g}\tfrom site A\t-1\t#0000ff\tnormal\t1\t0"
                               for g in ("gA", "gB")]


def test_text_dataset_rejects_a_tab_in_the_text(tmp_path):
    """A tab in the text would silently shift every field after it."""
    with pytest.raises(ItolError, match="tab"):
        write_text_dataset(str(tmp_path / "o.txt"), ["gA"], "a\tb")


################################################################################
# writes are atomic
################################################################################

def test_writers_leave_no_part_file_behind(tmp_path):
    out = tmp_path / "style.txt"
    write_style_dataset(str(out), ["gA"])
    assert out.exists()
    assert not os.path.exists(f"{out}.part")


################################################################################
# CLI flag vocabulary
################################################################################

def _subparsers():
    parser = build_parser()
    return next(a for a in parser._actions
                if isinstance(a, argparse._SubParsersAction)).choices


def test_no_short_flag_means_two_different_things():
    """
    Every short flag should carry one meaning across all subcommands. In bit these
    had drifted -- `-l` was `--label` under colorstrip and `--line-weight` under
    map, `-w` was `--width` and `--what-to-color` -- which silently does the wrong
    job when a command is copy-pasted between subcommands.
    """
    meanings = {}
    for sub in _subparsers().values():
        for action in sub._actions:
            for opt in action.option_strings:
                if opt.startswith("--") or opt in ("-h", "-v"):
                    continue
                meanings.setdefault(opt, set()).add(action.dest)

    collisions = {opt: sorted(d) for opt, d in meanings.items() if len(d) > 1}
    assert not collisions, f"short flags with more than one meaning: {collisions}"


def test_every_subcommand_takes_the_shared_arguments():
    for name, sub in _subparsers().items():
        dests = {a.dest for a in sub._actions}
        for shared in ("wanted_genomes", "input_dir", "source",
                       "output_file", "color", "dataset_label"):
            assert shared in dests, f"`{name}` is missing --{shared.replace('_', '-')}"


def test_check_args_rejects_a_non_numeric_line_weight():
    args = argparse.Namespace(func="style", line_weight="heavy",
                              wanted_genomes="x.txt", input_dir=None, source=None)
    with pytest.raises(ItolError, match="must be a number"):
        check_args(args)


def test_check_args_rejects_a_non_integer_strip_width():
    args = argparse.Namespace(func="colorstrip", strip_width="wide",
                              wanted_genomes="x.txt", input_dir=None, source=None)
    with pytest.raises(ItolError, match="must be an integer"):
        check_args(args)


################################################################################
# resolving a selection against a finished run
################################################################################

SUMMARY_ROWS = """genome_id\tinput\tsource\tlabel\ttaxid\tin_final_tree\treason_removed
GCF_000005845\tGCF_000005845.2\tncbi-accession (via -w)\tEscherichia_coli_GCF_000005845\t562\tYes\tNA
GCF_000009045\tGCF_000009045.1\tncbi-accession (via -w)\tBacillus_subtilis_GCF_000009045\t1423\tYes\tNA
my_MAG_01\tmy_MAG_01.fa\tnucleotide-fasta\tMy_MAG_01\tNA\tYes\tNA
my_MAG_02\tmy_MAG_02.fa\tnucleotide-fasta\tMy_MAG_02\tNA\tYes\tNA
my_MAG_03\tmy_MAG_03.fa\tnucleotide-fasta\tMy_MAG_03\tNA\tNo\ttoo few SCG hits
isolate_A\tisolate_A.gbk\tgenbank-file\tIsolate_A\tNA\tYes\tNA
"""


@pytest.fixture
def run_dir(tmp_path):
    d = tmp_path / "gtotree-output"
    d.mkdir()
    (d / "genomes-summary-info.tsv").write_text(SUMMARY_ROWS)
    return d


def test_load_summary_rejects_a_directory_that_isnt_a_run(tmp_path):
    from gtotree.utils.misc.itol import load_summary

    (tmp_path / "not-a-run").mkdir()
    with pytest.raises(ItolError, match="completed GToTree run"):
        load_summary(str(tmp_path / "not-a-run"))


def test_load_summary_rejects_a_missing_directory(tmp_path):
    from gtotree.utils.misc.itol import load_summary

    with pytest.raises(ItolError, match="not found"):
        load_summary(str(tmp_path / "nope"))


@pytest.mark.parametrize("identifier,expected", [
    ("my_MAG_01",       "My_MAG_01"),                       # genome_id
    ("My_MAG_01",       "My_MAG_01"),                       # label
    ("my_MAG_01.fa",    "My_MAG_01"),                       # original input
    ("GCF_000005845.2", "Escherichia_coli_GCF_000005845"),  # accession as given
    ("562",             "Escherichia_coli_GCF_000005845"),  # taxid
])
def test_wanted_genomes_resolve_from_any_identifier(identifier, expected, run_dir):
    """
    The list a user has on hand is whatever they kept -- accessions, the paths they
    passed in, IDs, or labels copied off the tree. All of them should land.
    """
    from gtotree.utils.misc.itol import resolve_wanted_genomes

    assert resolve_wanted_genomes(str(run_dir), [identifier]).labels == [expected]


def test_unmatched_genomes_are_a_hard_error(run_dir):
    """
    iToL silently ignores a row whose ID isn't in the tree, so a typo would
    otherwise show up as an annotation file that quietly does nothing.
    """
    from gtotree.utils.misc.itol import resolve_wanted_genomes

    with pytest.raises(ItolError, match="typo_genome"):
        resolve_wanted_genomes(str(run_dir), ["my_MAG_01", "typo_genome"])


def test_genomes_not_in_the_tree_are_dropped_and_reported(run_dir):
    from gtotree.utils.misc.itol import resolve_wanted_genomes

    selection = resolve_wanted_genomes(str(run_dir), ["my_MAG_01", "my_MAG_03"])
    assert selection.labels == ["My_MAG_01"]
    assert selection.dropped_not_in_tree == ["my_MAG_03"]


def test_selecting_only_genomes_absent_from_the_tree_raises(run_dir):
    from gtotree.utils.misc.itol import resolve_wanted_genomes

    with pytest.raises(ItolError, match="nothing to annotate"):
        resolve_wanted_genomes(str(run_dir), ["my_MAG_03"])


def test_ambiguous_identifiers_are_not_guessed_at(tmp_path):
    """
    A taxid shared by two genomes names neither of them, and resolving it either
    way would annotate the wrong leaf while looking correct.
    """
    from gtotree.utils.misc.itol import resolve_wanted_genomes

    d = tmp_path / "run"
    d.mkdir()
    (d / "genomes-summary-info.tsv").write_text(
        "genome_id\tinput\tsource\tlabel\ttaxid\tin_final_tree\treason_removed\n"
        "g1\tg1.fa\tnucleotide-fasta\tG1\t562\tYes\tNA\n"
        "g2\tg2.fa\tnucleotide-fasta\tG2\t562\tYes\tNA\n")

    with pytest.raises(ItolError, match="562"):
        resolve_wanted_genomes(str(d), ["562"])


################################################################################
# --source
################################################################################

def test_source_selects_without_any_id_list(run_dir):
    from gtotree.utils.misc.itol import select_by_source

    selection = select_by_source(str(run_dir), ["nucleotide-fasta"])
    assert selection.labels == ["My_MAG_01", "My_MAG_02"]
    assert selection.dropped_not_in_tree == ["my_MAG_03"]


def test_source_ncbi_accession_excludes_wanted_ref_tax_genomes(run_dir):
    """
    The whole point of selecting by source is usually separating the genomes you
    supplied from the backbone `-w` pulled in, so those must not be swept in by
    asking for `ncbi-accession`.
    """
    from gtotree.utils.misc.itol import select_by_source

    with pytest.raises(ItolError, match="No genomes"):
        select_by_source(str(run_dir), ["ncbi-accession"])

    selection = select_by_source(str(run_dir), ["wanted-ref-tax"])
    assert len(selection.labels) == 2


def test_multiple_sources_are_unioned(run_dir):
    from gtotree.utils.misc.itol import select_by_source

    selection = select_by_source(str(run_dir), ["nucleotide-fasta", "genbank-file"])
    assert set(selection.labels) == {"My_MAG_01", "My_MAG_02", "Isolate_A"}


################################################################################
# selection wiring
################################################################################

def _sel_args(**over):
    base = dict(wanted_genomes=None, input_dir=None, source=None)
    base.update(over)
    return argparse.Namespace(**base)


def test_check_args_requires_some_way_to_pick_genomes():
    with pytest.raises(ItolError, match="which genomes"):
        check_args(_sel_args(func="style", line_weight="3"))


def test_check_args_requires_a_run_dir_for_source():
    with pytest.raises(ItolError, match="`-i`"):
        check_args(_sel_args(func="style", line_weight="3",
                             source=["nucleotide-fasta"]))


def test_line_weight_loses_a_pointless_trailing_zero():
    """Keeps hand-made files byte-comparable with the auto-generated ones."""
    args = check_args(_sel_args(func="style", line_weight="3",
                                wanted_genomes="x.txt"))
    assert args.line_weight == 3


def test_source_and_wanted_genomes_intersect(run_dir, tmp_path):
    """`--source` sets the pool, `-g` narrows it."""
    from gtotree.cli.itol import select_targets

    wanted = tmp_path / "wanted.txt"
    wanted.write_text("my_MAG_01.fa\nisolate_A\n")

    targets, _dropped = select_targets(_sel_args(
        func="style", input_dir=str(run_dir), source=["nucleotide-fasta"],
        wanted_genomes=str(wanted)))

    assert targets == ["My_MAG_01"]


def test_wanted_genomes_without_a_run_dir_are_used_verbatim(tmp_path):
    """No run dir to consult means the list has to already be tree labels."""
    from gtotree.cli.itol import select_targets

    f = tmp_path / "labels.txt"
    f.write_text("Some_Label\nOther_Label\n")

    targets, dropped = select_targets(_sel_args(
        func="style", wanted_genomes=str(f)))

    assert targets == ["Some_Label", "Other_Label"]
    assert dropped == []


################################################################################
# automatic per-input-source files
################################################################################

def _summary(tmp_path, rows, name="run"):
    d = tmp_path / name
    d.mkdir()
    (d / "genomes-summary-info.tsv").write_text(
        "genome_id\tinput\tsource\tlabel\ttaxid\tin_final_tree\treason_removed\n"
        + rows)
    return str(d / "genomes-summary-info.tsv"), str(d / "itol-files")


MIXED = ("r1\tGCF_1\tncbi-accession (via -w)\tRef_1\tNA\tYes\tNA\n"
         "m1\tm1.fa\tnucleotide-fasta\tMy_MAG_01\tNA\tYes\tNA\n"
         "m2\tm2.fa\tnucleotide-fasta\tMy_MAG_02\tNA\tNo\ttoo few SCG hits\n"
         "i1\ti1.gbk\tgenbank-file\tIsolate_A\tNA\tYes\tNA\n")

ALL_FIVE = ("a1\tGCF_1\tncbi-accession\tAcc_1\tNA\tYes\tNA\n"
            "w1\tGCF_2\tncbi-accession (via -w)\tRef_1\tNA\tYes\tNA\n"
            "n1\tn1.fa\tnucleotide-fasta\tNuc_1\tNA\tYes\tNA\n"
            "p1\tp1.faa\tamino-acid-fasta\tAA_1\tNA\tYes\tNA\n"
            "g1\tg1.gbk\tgenbank-file\tGbk_1\tNA\tYes\tNA\n")


def test_a_file_is_written_for_every_source_present(tmp_path):
    from gtotree.utils.misc.itol import generate_input_source_itol_files

    summary, itol_dir = _summary(tmp_path, ALL_FIVE)
    written = generate_input_source_itol_files(summary, itol_dir)

    assert len(written) == 5
    assert sorted(os.listdir(itol_dir)) == [
        "amino-acid-fasta-itol.txt",
        "genbank-file-itol.txt",
        "ncbi-accession-itol.txt",
        "nucleotide-fasta-itol.txt",
        "wanted-ref-tax-itol.txt",
    ]


def test_wanted_ref_tax_file_is_slugged_not_named_after_its_label(tmp_path):
    """
    Its source label is "ncbi-accession (via -w)" -- spaces and parentheses make a
    miserable filename, and the slug matches the `--source` value that reproduces it.
    """
    from gtotree.utils.misc.itol import generate_input_source_itol_files

    summary, itol_dir = _summary(tmp_path, ALL_FIVE)
    generate_input_source_itol_files(summary, itol_dir)

    text = open(os.path.join(itol_dir, "wanted-ref-tax-itol.txt")).read()
    assert "DATASET_LABEL\tncbi-accession (via -w)\n" in text


def test_nucleotide_fasta_uses_the_same_blue_as_the_search_files(tmp_path):
    from gtotree.utils.misc.itol import (ITOL_COLOR,
                                         generate_input_source_itol_files)

    summary, itol_dir = _summary(tmp_path, ALL_FIVE)
    generate_input_source_itol_files(summary, itol_dir)

    text = open(os.path.join(itol_dir, "nucleotide-fasta-itol.txt")).read()
    assert f"COLOR\t{ITOL_COLOR}\n" in text


def test_every_source_file_gets_its_own_color(tmp_path):
    """They're meant to be loaded together, so shared colors would be useless."""
    from gtotree.utils.misc.itol import generate_input_source_itol_files

    summary, itol_dir = _summary(tmp_path, ALL_FIVE)
    generate_input_source_itol_files(summary, itol_dir)

    colors = set()
    for name in os.listdir(itol_dir):
        for line in open(os.path.join(itol_dir, name)):
            if line.startswith("COLOR\t"):
                colors.add(line.split("\t")[1].strip())
    assert len(colors) == 5


def test_source_files_skip_genomes_not_in_the_tree(tmp_path):
    from gtotree.utils.misc.itol import generate_input_source_itol_files

    summary, itol_dir = _summary(tmp_path, MIXED)
    generate_input_source_itol_files(summary, itol_dir)

    text = open(os.path.join(itol_dir, "nucleotide-fasta-itol.txt")).read()
    assert "My_MAG_01" in text
    assert "My_MAG_02" not in text


@pytest.mark.parametrize("rows,why", [
    ("a\ta.fa\tnucleotide-fasta\tA\tNA\tYes\tNA\n"
     "b\tb.fa\tnucleotide-fasta\tB\tNA\tYes\tNA\n",
     "one source -- every leaf would be the same color"),
    ("a\tGCF_1\tncbi-accession\tA\tNA\tYes\tNA\n"
     "b\tb.fa\tnucleotide-fasta\tB\tNA\tNo\ttoo few SCG hits\n",
     "the second source contributed nothing to the tree"),
])
def test_no_files_when_there_is_nothing_to_distinguish(rows, why, tmp_path):
    """These cases write nothing, and leave no empty directory behind either."""
    from gtotree.utils.misc.itol import generate_input_source_itol_files

    summary, itol_dir = _summary(tmp_path, rows)
    assert generate_input_source_itol_files(summary, itol_dir) == [], why
    assert not os.path.exists(itol_dir), why


def test_no_source_files_without_a_summary_table(tmp_path):
    from gtotree.utils.misc.itol import generate_input_source_itol_files

    assert generate_input_source_itol_files(
        str(tmp_path / "nope.tsv"), str(tmp_path / "out")) == []
