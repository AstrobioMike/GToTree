from gtotree.tests.smoke import main as smoke_main

def test_smoke_end_to_end(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert smoke_main([]) == 0
