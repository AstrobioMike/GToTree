"""
`gtt dl-ncbi-assemblies` -- taxon targeting, deduplication, and --dry-run.

The download engine itself is shared in spirit with bit's subcommand of the same name;
what's tested here is the GToTree-side wiring: that `-w` is repeatable, that `-a` and
`-w` combine into one deduplicated set, and that --dry-run reports without touching
anything.
"""

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest  # type: ignore

from gtotree.utils.ncbi.dl_ncbi_assemblies import (
    build_parser, resolve_targets, report_selection, dl_ncbi_assemblies,
    TaxonSelection, RunData)
from gtotree.utils.taxonomy.tax_select import CrossDomainTaxon, AmbiguousTaxon, TaxonNotFound

_MOD = "gtotree.utils.ncbi.dl_ncbi_assemblies"


def _args(**kwargs):
    defaults = dict(
        ncbi_accessions=None, wanted_ref_tax=None, source="gtdb",
        ncbi_section="refseq", derep_rank="auto", target_rank=None,
        target_domain=None, assembly_level=None, min_completeness=None,
        max_contamination=None, dry_run=False, format="fasta", jobs=10,
        output_dir=".",
    )
    defaults.update(kwargs)
    return SimpleNamespace(**defaults)


def _selection(canonical, accessions, rank="phylum", derep="family"):
    return SimpleNamespace(canonical=canonical, resolved_rank=rank,
                           effective_derep_rank=derep, accessions=list(accessions),
                           rows=[], warnings=[])


class TestParserFollowsGToTreeConventions:

    def test_w_is_the_taxon_flag_and_is_repeatable(self):
        args = build_parser().parse_args(["-w", "Nitrospirota", "-w", "Bacteria"])
        assert args.wanted_ref_tax == ["Nitrospirota", "Bacteria"]

    def test_a_is_the_accessions_flag(self):
        args = build_parser().parse_args(["-a", "accs.txt"])
        assert args.ncbi_accessions == "accs.txt"

    def test_source_defaults_to_gtdb_like_the_main_driver(self):
        assert build_parser().parse_args(["-w", "X"]).source == "gtdb"

    def test_derep_rank_defaults_to_auto(self):
        assert build_parser().parse_args(["-w", "X"]).derep_rank == "auto"

    def test_ncbi_section_defaults_to_both(self):
        """No section filter by default, so GenBank-only assemblies stay reachable."""
        assert build_parser().parse_args(["-w", "X"]).ncbi_section == "both"

    def test_fine_tuning_flags_are_hidden_until_detailed_help(self, capsys):
        build_parser().print_help()
        standard = capsys.readouterr().out
        build_parser(show_detailed=True).print_help()
        detailed = capsys.readouterr().out
        for flag in ("--min-completeness", "--max-contamination", "--assembly-level"):
            assert flag not in standard
            assert flag in detailed

    def test_hidden_flags_still_parse_and_keep_defaults(self):
        args = build_parser().parse_args(["-w", "X", "--min-completeness", "90"])
        assert args.min_completeness == 90.0
        assert build_parser().parse_args(["-w", "X"]).max_contamination is None


class TestResolveTargets:

    def test_accessions_file_alone(self, tmp_path):
        accs = tmp_path / "accs.txt"
        accs.write_text("GCA_001\nGCA_002\n")
        got, selections, (note, n_file) = resolve_targets(
            _args(ncbi_accessions=str(accs)))
        assert got == ["GCA_001", "GCA_002"]
        assert selections == [] and note is None and n_file == 2

    def test_repeated_w_dedupes_and_reports_the_overlap(self):
        first = _selection("Bacteria", ["GCF_1", "GCF_2"])
        second = _selection("Escherichia", ["GCF_2", "GCF_3"])
        with patch(f"{_MOD}.resolve_wanted_ref_tax_accessions",
                   side_effect=[(first.accessions, first), (second.accessions, second)]), \
             patch(f"{_MOD}.expand_wanted_ref_tax",
                   return_value=(["Bacteria", "Escherichia"], [])):
            got, selections, _ = resolve_targets(
                _args(wanted_ref_tax=["Bacteria", "Escherichia"]))

        assert got == ["GCF_1", "GCF_2", "GCF_3"]
        assert (selections[0].num_selected, selections[0].num_new) == (2, 2)
        assert (selections[1].num_selected, selections[1].num_new) == (2, 1)
        # the total is the union, not the sum of the per-taxon counts
        assert len(got) != sum(s.num_selected for s in selections)

    def test_w_dedupes_against_the_accessions_file(self, tmp_path):
        accs = tmp_path / "accs.txt"
        accs.write_text("GCF_1\n")
        sel = _selection("Bacteria", ["GCF_1", "GCF_2"])
        with patch(f"{_MOD}.resolve_wanted_ref_tax_accessions",
                   return_value=(sel.accessions, sel)), \
             patch(f"{_MOD}.expand_wanted_ref_tax", return_value=(["Bacteria"], [])):
            got, selections, _ = resolve_targets(
                _args(ncbi_accessions=str(accs), wanted_ref_tax=["Bacteria"]))
        assert got == ["GCF_1", "GCF_2"]
        assert selections[0].num_new == 1

    def test_ncbi_section_is_passed_to_the_resolver(self):
        sel = _selection("T", ["GCF_1"])
        with patch(f"{_MOD}.resolve_wanted_ref_tax_accessions",
                   return_value=(sel.accessions, sel)) as m, \
             patch(f"{_MOD}.expand_wanted_ref_tax", return_value=(["T"], [])):
            resolve_targets(_args(wanted_ref_tax=["T"], ncbi_section="genbank"))
        assert m.call_args.kwargs["ncbi_section"] == "genbank"

    def test_quality_floor_is_passed_to_the_resolver(self):
        sel = _selection("T", ["GCF_1"])
        with patch(f"{_MOD}.resolve_wanted_ref_tax_accessions",
                   return_value=(sel.accessions, sel)) as m, \
             patch(f"{_MOD}.expand_wanted_ref_tax", return_value=(["T"], [])):
            resolve_targets(_args(wanted_ref_tax=["T"], min_completeness=90.0,
                                  max_contamination=5.0))
        assert m.call_args.kwargs["min_completeness"] == 90.0
        assert m.call_args.kwargs["max_contamination"] == 5.0

    @pytest.mark.parametrize("err,expected", [
        (CrossDomainTaxon("Bacillus", ["Bacteria", "Eukaryota"]), "--target-domain"),
        (AmbiguousTaxon("X", ["order", "family"]), "--target-rank"),
        (TaxonNotFound("nope"), "doesn't seem to exist"),
    ])
    def test_resolution_errors_exit_with_gtotree_wording(self, err, expected, capsys):
        with patch(f"{_MOD}.expand_wanted_ref_tax", return_value=(["X"], [])), \
             patch(f"{_MOD}.resolve_wanted_ref_tax_accessions", side_effect=err):
            with pytest.raises(SystemExit):
                resolve_targets(_args(wanted_ref_tax=["X"]))
        assert expected in capsys.readouterr().out


class TestDryRun:

    def _run(self, tmp_path):
        first = _selection("Bacteria", ["GCF_1", "GCF_2"])
        second = _selection("Escherichia", ["GCF_2", "GCF_3"])
        a = _args(wanted_ref_tax=["Bacteria", "Escherichia"], dry_run=True,
                  output_dir=str(tmp_path / "out"))
        with patch(f"{_MOD}.resolve_wanted_ref_tax_accessions",
                   side_effect=[(first.accessions, first), (second.accessions, second)]), \
             patch(f"{_MOD}.expand_wanted_ref_tax",
                   return_value=(["Bacteria", "Escherichia"], [])), \
             patch(f"{_MOD}.get_ncbi_assembly_data"), \
             patch(f"{_MOD}.download_assemblies") as dl:
            dl_ncbi_assemblies(a)
        return dl, a

    def test_nothing_is_downloaded(self, tmp_path, capsys):
        dl, _ = self._run(tmp_path)
        capsys.readouterr()
        dl.assert_not_called()

    def test_no_output_dir_is_created(self, tmp_path, capsys):
        _, a = self._run(tmp_path)
        capsys.readouterr()
        assert not Path(a.output_dir).exists()

    def test_reports_per_taxon_counts_and_the_deduped_total(self, tmp_path, capsys):
        self._run(tmp_path)
        out = capsys.readouterr().out
        assert "Bacteria" in out and "Escherichia" in out
        assert "already counted" in out
        assert "3" in out          # union of {1,2} and {2,3}
        assert "dry run" in out.lower()


class TestReportSelection:

    def test_no_overlap_wording_when_there_is_no_overlap(self, capsys):
        sels = [TaxonSelection(taxon="A", canonical="A", resolved_rank="phylum",
                               effective_derep_rank="family", num_selected=5,
                               num_new=5, warnings=[])]
        report_selection(["a"] * 5, sels, (None, 0), _args())
        out = capsys.readouterr().out
        assert "5 genome(s) selected" in out
        assert "already counted" not in out

    def test_derep_off_is_labelled(self, capsys):
        sels = [TaxonSelection(taxon="A", canonical="A", resolved_rank="species",
                               effective_derep_rank=None, num_selected=2,
                               num_new=2, warnings=[])]
        report_selection(["a", "b"], sels, (None, 0), _args())
        assert "dereplication off" in capsys.readouterr().out

    def test_the_all_expansion_note_is_surfaced(self, capsys):
        sels = [TaxonSelection(taxon="Bacteria", canonical="Bacteria",
                               resolved_rank="domain", effective_derep_rank="class",
                               num_selected=1, num_new=1, warnings=[])]
        report_selection(["a"], sels,
                         ("`-w all` was expanded to: Bacteria, Archaea", 0), _args())
        assert "expanded to" in capsys.readouterr().out


class TestNotFoundMessagingBranchesOnInputMode:

    def test_accession_input_may_be_invalid(self):
        rd = RunData(from_taxon=False)
        assert "may be invalid" in rd.not_found_reason
        assert "assembly accessions?" in rd.none_found_hint

    def test_taxon_input_does_not_blame_the_user(self):
        rd = RunData(from_taxon=True)
        assert rd.not_found_reason == ""
        assert "invalid" not in rd.none_found_hint


class TestRegisteredInGtt:

    def test_subcommand_is_registered(self):
        from gtotree.cli.gtt import SUBCOMMAND_MAP
        assert SUBCOMMAND_MAP["dl-ncbi-assemblies"] == _MOD

    def test_listed_first_in_the_ncbi_gtdb_group(self):
        from gtotree.cli.gtt import PROGRAM_GROUPS
        group = next(g for g in PROGRAM_GROUPS if g["title"] == "NCBI/GTDB-related")
        assert group["programs"][0]["name"] == "dl-ncbi-assemblies"


class TestSelectionSkipsMetadataRows:
    """
    The downloader writes files, not a metadata TSV, and does no HMM auto-selection,
    so it asks for accessions only -- one less filtered read of the asset per taxon.
    The counts it reports still come from this same selection call, not from a cheaper
    distinct-group count (see
    test_domain_aware_derep.TestCheapCountPathWouldOverreport).
    """

    def _capture(self, **kwargs):
        sel = _selection("T", ["GCF_1"])
        with patch(f"{_MOD}.resolve_wanted_ref_tax_accessions",
                   return_value=(sel.accessions, sel)) as m, \
             patch(f"{_MOD}.expand_wanted_ref_tax", return_value=(["T"], [])):
            resolve_targets(_args(wanted_ref_tax=["T"], **kwargs))
        return m.call_args.kwargs

    def test_include_rows_is_false(self):
        assert self._capture()["include_rows"] is False

    def test_it_is_false_on_a_dry_run_too(self):
        assert self._capture(dry_run=True)["include_rows"] is False

    def test_dry_run_and_real_run_resolve_identically(self, tmp_path):
        sel = _selection("Bacteria", ["GCF_1", "GCF_2", "GCF_3"])

        def run(dry):
            a = _args(wanted_ref_tax=["Bacteria"], dry_run=dry,
                      output_dir=str(tmp_path / ("dry" if dry else "real")))
            with patch(f"{_MOD}.resolve_wanted_ref_tax_accessions",
                       return_value=(sel.accessions, sel)), \
                 patch(f"{_MOD}.expand_wanted_ref_tax",
                       return_value=(["Bacteria"], [])):
                return resolve_targets(a)

        dry_accs, dry_sels, _ = run(True)
        real_accs, real_sels, _ = run(False)
        assert dry_accs == real_accs
        assert dry_sels[0].num_selected == real_sels[0].num_selected


################################################################################
# download_one: built on urllib, not requests.
#
# GToTree has no `requests` dependency (nothing else in the package imports it, and
# it isn't in the conda recipe), so this subcommand uses urllib like the rest of the
# codebase. The retry/backoff POLICY is imported from processing_genomes rather than
# duplicated.
################################################################################

import io
import urllib.error
from unittest.mock import MagicMock

from gtotree.utils.ncbi import dl_ncbi_assemblies as dlmod


def _resp(body=b"data", content_type="application/octet-stream"):
    """A stand-in for the context-manager object urlopen returns."""
    resp = MagicMock()
    resp.headers = {"Content-Type": content_type}
    resp.__enter__.return_value = resp
    resp.__exit__.return_value = False
    resp.read = io.BytesIO(body).read
    return resp


def _http_error(code, headers=None):
    return urllib.error.HTTPError("http://x", code, "err", headers or {}, None)


class TestDownloaderUsesUrllibNotRequests:

    def test_the_module_does_not_import_requests(self):
        import inspect
        src = inspect.getsource(dlmod)
        assert "import requests" not in src
        assert "urllib.request" in src

    def test_it_reuses_the_shared_backoff_policy(self):
        """The subtle part (throttle split, sawtooth, ceilings) is shared, not copied."""
        import inspect
        src = inspect.getsource(dlmod)
        assert "from gtotree.utils.misc.processing_genomes import" in src
        assert "_sleep_backoff" in src


class TestDownloadOne:

    def test_a_good_response_is_written_atomically(self, tmp_path, monkeypatch):
        dest = tmp_path / "GCF_1.fasta.gz"
        monkeypatch.setattr(dlmod.urllib.request, "urlopen",
                            lambda *a, **k: _resp(b"payload"))
        result = dlmod.download_one("http://x/f.gz", str(dest))
        assert result == (str(dest), None, "downloaded")
        assert dest.read_bytes() == b"payload"
        assert not (tmp_path / "GCF_1.fasta.gz.tmp").exists()

    def test_an_existing_nonempty_file_is_skipped(self, tmp_path, monkeypatch):
        dest = tmp_path / "GCF_1.fasta.gz"
        dest.write_bytes(b"already here")
        def boom(*a, **k):
            raise AssertionError("should not have hit the network")
        monkeypatch.setattr(dlmod.urllib.request, "urlopen", boom)
        assert dlmod.download_one("http://x/f.gz", str(dest))[2] == "skipped"

    def test_404_is_permanent_and_not_retried(self, tmp_path, monkeypatch):
        calls = []
        def urlopen(*a, **k):
            calls.append(1)
            raise _http_error(404)
        monkeypatch.setattr(dlmod.urllib.request, "urlopen", urlopen)
        dest = tmp_path / "GCF_1.fasta.gz"
        _d, err, status = dlmod.download_one("http://x/f.gz", str(dest))
        assert status == "failed"
        assert "404" in err
        assert len(calls) == 1, "a 404 must not be retried"

    def test_a_500_is_transient_and_retried(self, tmp_path, monkeypatch):
        calls = []
        def urlopen(*a, **k):
            calls.append(1)
            raise _http_error(500)
        monkeypatch.setattr(dlmod.urllib.request, "urlopen", urlopen)
        monkeypatch.setattr(dlmod, "_sleep_backoff", lambda *a, **k: None)
        dest = tmp_path / "GCF_1.fasta.gz"
        _d, _err, status = dlmod.download_one("http://x/f.gz", str(dest), retries=3)
        assert status == "failed_transient"
        assert len(calls) == 3

    def test_a_429_is_treated_as_throttling(self, tmp_path, monkeypatch):
        seen = {}
        monkeypatch.setattr(dlmod.urllib.request, "urlopen",
                            lambda *a, **k: (_ for _ in ()).throw(_http_error(429)))
        monkeypatch.setattr(dlmod, "_sleep_backoff",
                            lambda attempt, err=None, throttled=False: seen.update(
                                throttled=throttled))
        dest = tmp_path / "GCF_1.fasta.gz"
        dlmod.download_one("http://x/f.gz", str(dest), retries=2)
        assert seen.get("throttled") is True

    def test_an_html_error_page_is_not_written_as_the_file(self, tmp_path, monkeypatch):
        monkeypatch.setattr(dlmod.urllib.request, "urlopen",
                            lambda *a, **k: _resp(b"<html>oops</html>", "text/html"))
        monkeypatch.setattr(dlmod, "_sleep_backoff", lambda *a, **k: None)
        dest = tmp_path / "GCF_1.fasta.gz"
        _d, err, status = dlmod.download_one("http://x/f.gz", str(dest), retries=2)
        assert status == "failed_transient"
        assert "error page" in err
        assert not dest.exists(), "an error page must never land as the real file"

    def test_an_empty_body_is_not_left_behind(self, tmp_path, monkeypatch):
        monkeypatch.setattr(dlmod.urllib.request, "urlopen",
                            lambda *a, **k: _resp(b""))
        monkeypatch.setattr(dlmod, "_sleep_backoff", lambda *a, **k: None)
        dest = tmp_path / "GCF_1.fasta.gz"
        _d, err, status = dlmod.download_one("http://x/f.gz", str(dest), retries=2)
        assert status == "failed_transient"
        assert "empty" in err
        assert not dest.exists()
        assert not (tmp_path / "GCF_1.fasta.gz.tmp").exists()

    def test_a_transport_error_is_transient(self, tmp_path, monkeypatch):
        monkeypatch.setattr(dlmod.urllib.request, "urlopen",
                            lambda *a, **k: (_ for _ in ()).throw(
                                urllib.error.URLError("connection reset")))
        monkeypatch.setattr(dlmod, "_sleep_backoff", lambda *a, **k: None)
        dest = tmp_path / "GCF_1.fasta.gz"
        assert dlmod.download_one("http://x/f.gz", str(dest), retries=2)[2] == \
            "failed_transient"
