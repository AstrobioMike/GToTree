from gtotree.tests.smoke import run_smoke_test

def test_smoke_end_to_end(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert run_smoke_test([]) == 0
