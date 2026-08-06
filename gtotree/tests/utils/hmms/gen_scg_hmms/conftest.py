"""
Fixtures for the `gtt gen-scg-hmms` tests.
"""

import pytest  # type: ignore


@pytest.fixture(autouse=True)
def reference_data_fetches(monkeypatch):
    """
    Stub the managed-dataset fetches, and record which ones were asked for.

    Autouse because `phase_resolve_genomes` now ensures the NCBI/GTDB tables are on
    disk before anything reads them, so any test that resolves accessions would
    otherwise try to download a multi-hundred-MB asset. That's the kind of thing you
    want impossible in a test package rather than remembered case by case.

    Requestable by name too, so a test can assert on what got fetched. Doing it with one
    fixture rather than an autouse stub plus a separate tracker avoids depending on
    fixture ordering to decide which patch wins.

    Patched on the defining modules because the helper imports them inside its body,
    which resolves the attribute at call time.
    """
    import gtotree.utils.ncbi.get_ncbi_assembly_data as ncbi_mod
    import gtotree.utils.gtdb.get_gtdb_data as gtdb_mod

    calls = []
    monkeypatch.setattr(ncbi_mod, "get_ncbi_assembly_data",
                        lambda *a, **kw: calls.append("ncbi"))
    monkeypatch.setattr(gtdb_mod, "get_gtdb_data",
                        lambda *a, **kw: calls.append("gtdb"))
    return calls
