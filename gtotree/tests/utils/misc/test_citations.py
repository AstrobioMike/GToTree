"""
Unit tests for gtotree/utils/misc/citations.py.

citations.txt is what people paste into a manuscript's methods section, so the two
things that matter are that the tools a run actually used are the ones listed, and that
the versions quoted are the versions that ran. Both are easy to get subtly wrong and
neither fails loudly: a missing citation or a wrong version number just sits there
looking plausible.

The version strings below are real `--version` output from the pinned tools in
conda-recipe/meta.yaml. They're the reason this file exists -- `_extract_version` is a
regex over banner text, which is exactly the kind of thing that silently starts
returning "" when a tool reformats its banner.
"""

import subprocess

import pytest  # type: ignore

from gtotree.utils.misc import citations
from gtotree.utils.misc.citations import (CitationsInfo, _extract_version,
                                          _tool_output, generate_citations_info,
                                          get_hmmer_label)
from gtotree.utils.misc.general import ToolsUsed


#: real banner output, keyed by the getter that parses it
TOOL_BANNERS = {
    "muscle": "muscle 5.1.linux64 []",
    "trimal": "trimAl v1.4.rev15 build[2013-12-17]",
    "prodigal": "\nProdigal V2.6.3: February, 2016\n",
    "fasttree": "Detailed usage for FastTree 2.1.11 Double precision (No SSE3)",
    "veryfasttree": "VeryFastTree Version 4.0.3 (Built with cmake)",
    "iqtree": "IQ-TREE multicore version 2.2.0 COVID-edition for Linux 64-bit",
    "kofamscan": "exec_annotation 1.3.0",
}


class _RunData:
    def __init__(self, output_dir, **used):
        self.output_dir = str(output_dir)
        self.tools_used = ToolsUsed(**used)


# ---------------------------------------------------------------------------
# _extract_version
# ---------------------------------------------------------------------------

class TestExtractVersion:

    @pytest.mark.parametrize("tool,expected", [
        ("muscle", "5.1.linux64"),
        ("trimal", "1.4.rev15"),
        ("prodigal", "2.6.3"),
        ("fasttree", "2.1.11"),
        ("veryfasttree", "4.0.3"),
        ("iqtree", "2.2.0"),
        ("kofamscan", "1.3.0"),
    ])
    def test_a_real_banner_yields_the_version(self, tool, expected):
        assert _extract_version(TOOL_BANNERS[tool]) == expected

    def test_a_marker_narrows_to_the_right_line(self):
        # a banner that mentions another version before its own
        text = "built against zlib 1.2.11\ntrimAl v1.4.rev15 build[2013-12-17]"
        assert _extract_version(text, marker="trim") == "1.4.rev15"

    def test_the_marker_is_case_insensitive(self):
        assert _extract_version(TOOL_BANNERS["prodigal"], marker="PRODIGAL") == "2.6.3"

    def test_a_marker_that_matches_nothing_falls_back_to_every_line(self):
        """
        Deliberate: a tool that reformats its banner and drops the marker word should
        still get a version out, rather than silently reporting none.
        """
        assert _extract_version(TOOL_BANNERS["iqtree"], marker="nonsense") == "2.2.0"

    def test_the_first_version_shaped_token_wins(self):
        assert _extract_version("tool 1.2.3\nlibrary 9.9.9") == "1.2.3"

    def test_a_bare_integer_is_not_a_version(self):
        # the pattern needs at least one dot-separated part, or every year and count
        # in a banner would match
        assert _extract_version("some tool build 2016") == ""

    @pytest.mark.parametrize("text", ["", "no version here", "\n\n"])
    def test_text_without_a_version_yields_empty(self, text):
        assert _extract_version(text) == ""


# ---------------------------------------------------------------------------
# _tool_output and the per-tool getters
# ---------------------------------------------------------------------------

class TestToolOutput:

    def test_stdout_and_stderr_are_both_searched(self, monkeypatch):
        # tools are split on which stream they print their banner to
        def _fake(cmd, capture_output, text):
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="tool 1.2.3")
        monkeypatch.setattr(citations.subprocess, "run", _fake)

        assert "1.2.3" in _tool_output(["whatever"])

    def test_a_missing_binary_is_empty_rather_than_an_exception(self, monkeypatch):
        """
        citations.txt is written at the very end of a successful run. A tool that
        isn't on PATH must not take the run down at that point.
        """
        def _boom(*a, **k):
            raise FileNotFoundError("no such binary")
        monkeypatch.setattr(citations.subprocess, "run", _boom)

        assert _tool_output(["not-a-real-tool"]) == ""


class TestVersionGetters:

    GETTERS = {
        "muscle": ("get_muscle_version", "5.1.linux64"),
        "trimal": ("get_trimal_version", "v1.4.rev15"),
        "prodigal": ("get_prodigal_version", "v2.6.3"),
        "fasttree": ("get_fasttree_version", "2.1.11"),
        "veryfasttree": ("get_veryfasttree_version", "4.0.3"),
        "iqtree": ("get_iqtree_version", "2.2.0"),
        "kofamscan": ("get_kofamscan_version", "1.3.0"),
    }

    @pytest.mark.parametrize("tool", sorted(GETTERS))
    def test_each_getter_parses_its_tools_banner(self, tool, monkeypatch):
        name, expected = self.GETTERS[tool]
        monkeypatch.setattr(citations, "_tool_output",
                            lambda cmd: TOOL_BANNERS[tool])
        assert getattr(citations, name)() == expected

    @pytest.mark.parametrize("tool", sorted(GETTERS))
    def test_each_getter_is_empty_when_its_tool_is_absent(self, tool, monkeypatch):
        # and specifically not the bare "v" that a naive f-string prefix would give
        name, _ = self.GETTERS[tool]
        monkeypatch.setattr(citations, "_tool_output", lambda cmd: "")
        assert getattr(citations, name)() == ""


class TestHmmerLabel:

    def test_the_hmmer_version_is_preferred(self, monkeypatch):
        monkeypatch.setattr(citations, "get_hmmer_version", lambda: "3.4")
        assert get_hmmer_label() == "HMMER v3.4"

    def test_it_falls_back_to_the_pyhmmer_version(self, monkeypatch):
        import pyhmmer  # type: ignore
        monkeypatch.setattr(citations, "get_hmmer_version", lambda: None)
        assert get_hmmer_label() == f"HMMER (via pyhmmer v{pyhmmer.__version__})"


# ---------------------------------------------------------------------------
# generate_citations_info
# ---------------------------------------------------------------------------

class TestGenerateCitationsInfo:
    """
    The part that would be embarrassing to get wrong: citing a tool that never ran, or
    failing to cite one that did.
    """

    @pytest.fixture(autouse=True)
    def stub_versions(self, monkeypatch):
        # no external binaries in the test environment; the parsing is covered above
        monkeypatch.setattr(citations, "_tool_output", lambda cmd: "tool 9.9.9")
        # and don't depend on whether the package is pip-installed in this environment
        monkeypatch.setattr(citations, "version", lambda name: "2.0.0")

    def test_the_gtotree_version_comes_from_the_installed_distribution(self,
                                                                      tmp_path,
                                                                      monkeypatch):
        """
        The lookup is by distribution name, which importlib.metadata normalizes -- so
        "GToTree" resolves against the "gtotree" in pyproject.toml. Pinning it here
        because a rename in pyproject would otherwise surface as a PackageNotFoundError
        at the very end of an otherwise successful run.
        """
        seen = []
        monkeypatch.setattr(citations, "version",
                            lambda name: seen.append(name) or "2.0.0")

        generate_citations_info(_RunData(tmp_path))

        assert seen == ["GToTree"]
        assert "GToTree v2.0.0" in (tmp_path / "citations.txt").read_text()

    def _generate(self, tmp_path, **used):
        generate_citations_info(_RunData(tmp_path, **used))
        return (tmp_path / "citations.txt").read_text()

    def test_the_always_on_tools_are_always_cited(self, tmp_path):
        text = self._generate(tmp_path)

        # GToTree, HMMER, MUSCLE and TrimAl run on every single tree
        assert CitationsInfo.gtotree in text
        assert CitationsInfo.hmmer in text
        assert CitationsInfo.muscle in text
        assert CitationsInfo.trimal in text

    @pytest.mark.parametrize("flag,attr", [
        ("prodigal_used", "prodigal"),
        ("gtdb_used", "gtdb"),
        ("fasttree_used", "fasttree"),
        ("veryfasttree_used", "veryfasttree"),
        ("iqtree_used", "iqtree"),
        ("kofamscan_used", "kofamscan"),
        ("pfam_db_used", "pfam"),
        ("universal_SCGs_used", "universal_SCG_set"),
    ])
    def test_an_optional_tool_is_cited_only_when_it_was_used(self, tmp_path, flag,
                                                             attr):
        citation = getattr(CitationsInfo, attr)

        assert citation not in self._generate(tmp_path)
        assert citation in self._generate(tmp_path, **{flag: True})

    def test_only_the_tree_program_that_ran_is_cited(self, tmp_path):
        # the three are mutually exclusive in a run; citing two would be wrong
        text = self._generate(tmp_path, fasttree_used=True)

        assert CitationsInfo.fasttree in text
        assert CitationsInfo.veryfasttree not in text
        assert CitationsInfo.iqtree not in text

    def test_the_file_lands_in_the_output_dir(self, tmp_path):
        generate_citations_info(_RunData(tmp_path))
        assert (tmp_path / "citations.txt").is_file()

    def test_rerunning_overwrites_rather_than_appends(self, tmp_path):
        # resume re-runs this stage; a doubled citations file would be a visible bug
        first = self._generate(tmp_path, prodigal_used=True)
        second = self._generate(tmp_path, prodigal_used=True)

        assert first == second
        assert second.count(CitationsInfo.prodigal) == 1

    def test_a_version_that_cannot_be_read_still_writes_the_citation(self, tmp_path,
                                                                     monkeypatch):
        """
        A tool whose banner we can't parse should cost us the version number, not the
        reference -- the citation is the part that matters for a methods section.
        """
        monkeypatch.setattr(citations, "_tool_output", lambda cmd: "")

        text = self._generate(tmp_path, prodigal_used=True)

        assert CitationsInfo.prodigal in text
        assert CitationsInfo.muscle in text
