"""
Unit tests for gtotree/utils/misc/data_locations.py.

This is the module behind `gtt data locations check` and `gtt data locations set`, so
it runs on essentially every fresh install and it WRITES into the user's conda env
(etc/conda/activate.d/gtotree.sh), or into ~/.config/gtotree/ when they can't write
there. A bug here doesn't fail a run, it mangles someone's environment -- which is why
the rewrite logic gets the most attention below.

The rewrite contract, which is easy to break by accident:

  * lines for variables we're about to set are DROPPED, then re-appended, so running
    `set` twice doesn't leave two exports for the same variable
  * lines for anything else (LC_ALL, LANG, and whatever a user added) are KEPT
  * when falling back to ~/.config/gtotree/, the `if ... . ~/.config/gtotree/gtotree.sh
    ... fi` block from build.sh must NOT be copied across, or the file sources itself

Nothing here touches a real conda env -- CONDA_PREFIX is pointed at a tmp_path.
"""

import os

import pytest  # type: ignore

from gtotree.utils.misc import data_locations as dl


ACTIVATE_REL = "etc/conda/activate.d/gtotree.sh"

# what conda-recipe/build.sh actually lays down, trimmed to two variables
BUILD_SH_SCRIPT = """\
export NCBI_ASSEMBLY_DATA_DIR=${CONDA_PREFIX}/share/gtotree/ncbi_assembly_summaries/
export GTDB_DIR=${CONDA_PREFIX}/share/gtotree/gtdb_tax_info/
export LC_ALL="en_US.UTF-8"
export LANG="en_US.UTF-8"

if [ -f ~/.config/gtotree/gtotree.sh ] && [ ! -w ${CONDA_PREFIX}/etc/conda/activate.d/gtotree.sh ]; then
    . ~/.config/gtotree/gtotree.sh
fi
"""


@pytest.fixture
def conda_env(tmp_path, monkeypatch):
    """A fake conda prefix with build.sh's activate script already in place."""
    script = tmp_path / ACTIVATE_REL
    script.parent.mkdir(parents=True)
    script.write_text(BUILD_SH_SCRIPT)
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
    return script


@pytest.fixture
def fake_home(tmp_path, monkeypatch):
    """Point ~ at a tmp dir so the not-writable fallback can't touch a real home."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(os.path, "expanduser", lambda p: p.replace("~", str(home), 1))
    return home


def _unwritable(monkeypatch, *paths):
    """
    Make `os.access(path, W_OK)` report False for `paths`, leaving everything else
    alone. Done by patching rather than chmod because root ignores the permission
    bits entirely, which would silently turn these into no-op tests in a container.
    """
    targets = {str(p) for p in paths}
    real_access = os.access

    def fake_access(path, mode, **kwargs):
        if str(path) in targets and mode == os.W_OK:
            return False
        return real_access(path, mode, **kwargs)

    monkeypatch.setattr(dl.os, "access", fake_access)


def _exports(text):
    """{variable: value} for every `export VAR=value` line in a script."""
    out = {}
    for line in text.splitlines():
        if line.startswith("export "):
            var, _, value = line[len("export "):].partition("=")
            out[var] = value
    return out


# ---------------------------------------------------------------------------
# check_location_var_is_set_and_writable
# ---------------------------------------------------------------------------

class TestCheckLocationVar:

    def test_a_set_writable_path_comes_back_writable(self, tmp_path, monkeypatch):
        monkeypatch.setenv("GTDB_DIR", str(tmp_path))
        assert dl.check_location_var_is_set_and_writable("GTDB_DIR") == \
            (str(tmp_path), True)

    def test_a_set_but_unwritable_path_is_reported_not_fatal(self, tmp_path,
                                                             monkeypatch):
        """
        An unwritable location is a warning, not an exit -- `locations check` has to
        get through every variable to print the table.
        """
        locked = tmp_path / "locked"
        locked.mkdir()
        _unwritable(monkeypatch, locked)
        monkeypatch.setenv("GTDB_DIR", str(locked))

        assert dl.check_location_var_is_set_and_writable("GTDB_DIR") == \
            (str(locked), False)

    def test_an_unset_variable_exits(self, monkeypatch, capsys):
        monkeypatch.delenv("GTDB_DIR", raising=False)
        with pytest.raises(SystemExit) as excinfo:
            dl.check_location_var_is_set_and_writable("GTDB_DIR")
        assert excinfo.value.code == 1
        assert "does not seem to be set" in capsys.readouterr().out

    def test_an_empty_variable_counts_as_unset(self, monkeypatch):
        monkeypatch.setenv("GTDB_DIR", "")
        with pytest.raises(SystemExit):
            dl.check_location_var_is_set_and_writable("GTDB_DIR")


class TestCheckAndReportEnvVariables:

    def test_every_variable_is_listed(self, tmp_path, monkeypatch, capsys):
        for variable in dl.ENV_VARIABLES:
            monkeypatch.setenv(variable, str(tmp_path))

        dl.check_and_report_env_variables()

        out = capsys.readouterr().out
        for variable in dl.ENV_VARIABLES:
            assert variable in out

    def test_an_unwritable_location_is_called_out(self, tmp_path, monkeypatch,
                                                  capsys):
        locked = tmp_path / "locked"
        locked.mkdir()
        _unwritable(monkeypatch, locked)
        for variable in dl.ENV_VARIABLES:
            monkeypatch.setenv(variable, str(tmp_path))
        monkeypatch.setenv("Pfam_data_dir", str(locked))

        dl.check_and_report_env_variables()

        out = capsys.readouterr().out
        assert "'Pfam_data_dir' variable is not writable" in out


class TestGetVariablePath:

    def test_returns_the_path_when_set(self, monkeypatch):
        monkeypatch.setenv("KO_data_dir", "/somewhere")
        assert dl.get_variable_path("KO_data_dir") == "/somewhere"

    def test_returns_false_when_unset(self, monkeypatch):
        monkeypatch.delenv("KO_data_dir", raising=False)
        assert dl.get_variable_path("KO_data_dir") is False


# ---------------------------------------------------------------------------
# set_variable_path -- the interactive prompt
# ---------------------------------------------------------------------------

class TestSetVariablePath:

    def test_declining_to_change_keeps_the_current_path(self, monkeypatch):
        monkeypatch.setattr("builtins.input", lambda _: "n")
        assert dl.set_variable_path("GTDB_DIR", "/current/path") == "/current/path"

    def test_a_junk_answer_is_re_asked(self, monkeypatch, capsys):
        answers = iter(["maybe", "n"])
        monkeypatch.setattr("builtins.input", lambda _: next(answers))

        assert dl.set_variable_path("GTDB_DIR", "/current/path") == "/current/path"
        assert "Must respond with 'y' or 'n'" in capsys.readouterr().out

    def test_an_unset_variable_goes_straight_to_the_prompt(self, tmp_path,
                                                          monkeypatch):
        wanted = tmp_path / "new-spot"
        monkeypatch.setattr("builtins.input", lambda _: str(wanted))

        result = dl.set_variable_path("GTDB_DIR", False)

        assert result == str(wanted) + os.sep
        assert wanted.is_dir()  # a path that doesn't exist yet is created

    def test_agreeing_to_change_then_giving_a_path(self, tmp_path, monkeypatch):
        wanted = tmp_path / "new-spot"
        answers = iter(["y", str(wanted)])
        monkeypatch.setattr("builtins.input", lambda _: next(answers))

        assert dl.set_variable_path("GTDB_DIR", "/current/path") == \
            str(wanted) + os.sep

    def test_an_unwritable_path_is_re_asked_rather_than_accepted(self, tmp_path,
                                                                 monkeypatch, capsys):
        locked = tmp_path / "locked"
        locked.mkdir()
        # os.path.join(path, "") appends a separator before the access check
        _unwritable(monkeypatch, str(locked) + os.sep)
        good = tmp_path / "good"
        answers = iter([str(locked), str(good)])
        monkeypatch.setattr("builtins.input", lambda _: next(answers))

        result = dl.set_variable_path("GTDB_DIR", False)

        assert result == str(good) + os.sep
        assert "not writable for you" in capsys.readouterr().out

    def test_a_relative_path_is_re_asked_rather_than_accepted(self, tmp_path,
                                                              monkeypatch, capsys):
        monkeypatch.chdir(tmp_path)
        good = tmp_path / "good"
        answers = iter(["relative-spot", str(good)])
        monkeypatch.setattr("builtins.input", lambda _: next(answers))

        result = dl.set_variable_path("GTDB_DIR", False)

        assert result == str(good) + os.sep
        assert "absolute path" in capsys.readouterr().out


class TestSetEnvVariables:

    def test_collects_a_path_for_every_managed_variable(self, tmp_path, monkeypatch):
        for variable in dl.ENV_VARIABLES:
            monkeypatch.delenv(variable, raising=False)
        monkeypatch.setattr("builtins.input", lambda _: str(tmp_path / "spot"))

        paths = dl.set_env_variables()

        assert sorted(paths) == sorted(dl.ENV_VARIABLES)


# ---------------------------------------------------------------------------
# modify_conda_activate_startup_script -- the part that writes to disk
# ---------------------------------------------------------------------------

class TestRewritingTheWritableCondaScript:

    def test_managed_variables_are_replaced_not_duplicated(self, conda_env):
        dl.modify_conda_activate_startup_script({"GTDB_DIR": "/new/gtdb"})

        text = conda_env.read_text()
        assert text.count("export GTDB_DIR=") == 1
        assert _exports(text)["GTDB_DIR"] == "/new/gtdb"

    def test_unmanaged_lines_are_preserved(self, conda_env):
        dl.modify_conda_activate_startup_script({"GTDB_DIR": "/new/gtdb"})

        exports = _exports(conda_env.read_text())
        assert exports["LC_ALL"] == '"en_US.UTF-8"'
        assert exports["LANG"] == '"en_US.UTF-8"'

    def test_variables_not_being_set_are_left_alone(self, conda_env):
        dl.modify_conda_activate_startup_script({"GTDB_DIR": "/new/gtdb"})

        exports = _exports(conda_env.read_text())
        assert exports["NCBI_ASSEMBLY_DATA_DIR"] == \
            "${CONDA_PREFIX}/share/gtotree/ncbi_assembly_summaries/"

    def test_the_home_sourcing_block_survives_in_the_conda_script(self, conda_env):
        # it's only skipped when writing the home-location copy, not here
        dl.modify_conda_activate_startup_script({"GTDB_DIR": "/new/gtdb"})
        assert ". ~/.config/gtotree/gtotree.sh" in conda_env.read_text()

    def test_running_twice_is_idempotent(self, conda_env):
        dl.modify_conda_activate_startup_script({"GTDB_DIR": "/new/gtdb"})
        first = conda_env.read_text()
        dl.modify_conda_activate_startup_script({"GTDB_DIR": "/new/gtdb"})

        assert conda_env.read_text() == first

    def test_no_tmp_file_is_left_behind(self, conda_env):
        dl.modify_conda_activate_startup_script({"GTDB_DIR": "/new/gtdb"})
        assert not (conda_env.parent / "gtotree.sh.tmp").exists()


class TestFallingBackToTheHomeLocation:

    @pytest.fixture
    def unwritable_conda(self, conda_env, monkeypatch):
        _unwritable(monkeypatch, conda_env)
        return conda_env

    def test_the_conda_script_is_not_touched(self, unwritable_conda, fake_home):
        before = unwritable_conda.read_text()
        dl.modify_conda_activate_startup_script({"GTDB_DIR": "/new/gtdb"})
        assert unwritable_conda.read_text() == before

    def test_the_settings_land_in_the_home_config(self, unwritable_conda, fake_home):
        dl.modify_conda_activate_startup_script({"GTDB_DIR": "/new/gtdb"})

        written = fake_home / ".config/gtotree/gtotree.sh"
        assert written.exists()
        assert _exports(written.read_text())["GTDB_DIR"] == "/new/gtdb"

    def test_the_self_sourcing_block_is_not_copied_across(self, unwritable_conda,
                                                          fake_home):
        """
        The home copy is what the conda script sources. Carrying the sourcing block
        into it would make it source itself.
        """
        dl.modify_conda_activate_startup_script({"GTDB_DIR": "/new/gtdb"})

        text = (fake_home / ".config/gtotree/gtotree.sh").read_text()
        assert ". ~/.config/gtotree/gtotree.sh" not in text
        assert "fi" not in text


class TestMissingOrAbsentEnvironment:
    """
    The two first-run shapes that used to raise a bare traceback out of `gtt data
    locations set` rather than saying anything useful.
    """

    def test_a_conda_env_with_no_startup_script_yet_is_written_fresh(
            self, tmp_path, monkeypatch, fake_home):
        # os.access() is False for a missing file, so this takes the home-fallback
        # path -- it must not then try to read the script that isn't there
        (tmp_path / "etc/conda/activate.d").mkdir(parents=True)
        monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))

        dl.modify_conda_activate_startup_script({"GTDB_DIR": "/new/gtdb"})

        written = fake_home / ".config/gtotree/gtotree.sh"
        assert _exports(written.read_text()) == {"GTDB_DIR": "/new/gtdb"}

    def test_outside_a_conda_env_it_explains_itself_and_exits(self, monkeypatch,
                                                              capsys):
        monkeypatch.delenv("CONDA_PREFIX", raising=False)

        with pytest.raises(SystemExit) as excinfo:
            dl.modify_conda_activate_startup_script({"GTDB_DIR": "/new/gtdb"})

        assert excinfo.value.code == 1
        out = capsys.readouterr().out
        assert "doesn't look like a conda environment" in out
        # the variables are named so a non-conda user can set them by hand
        assert "GTDB_DIR" in out

    def test_nowhere_writable_at_all_points_at_the_issue_tracker(
            self, conda_env, tmp_path, monkeypatch, capsys):
        _unwritable(monkeypatch, conda_env)
        monkeypatch.setattr(dl.os, "makedirs",
                            lambda *a, **k: (_ for _ in ()).throw(OSError("nope")))

        with pytest.raises(SystemExit) as excinfo:
            dl.modify_conda_activate_startup_script({"GTDB_DIR": "/new/gtdb"})

        assert excinfo.value.code == 1
        assert "GToTree/issues" in capsys.readouterr().out


class TestNotifyToReactivateConda:

    def test_names_the_env_to_reactivate(self, monkeypatch, capsys):
        monkeypatch.setenv("CONDA_DEFAULT_ENV", "gtotree-dev")
        dl.notify_to_reactivate_conda()
        assert "conda deactivate && conda activate gtotree-dev" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# ensure_reference_data -- which assets a run actually pulls
# ---------------------------------------------------------------------------

class TestEnsureReferenceData:

    @pytest.fixture
    def fetched(self, monkeypatch):
        calls = []
        monkeypatch.setattr(
            "gtotree.utils.ncbi.get_ncbi_assembly_data.get_ncbi_assembly_data",
            lambda *a, **k: calls.append("ncbi"))
        monkeypatch.setattr("gtotree.utils.gtdb.get_gtdb_data.get_gtdb_data",
                            lambda *a, **k: calls.append("gtdb"))
        return calls

    def test_nothing_is_fetched_for_a_run_that_needs_neither(self, fetched):
        dl.ensure_reference_data()
        assert fetched == []

    def test_ncbi_accessions_pull_the_ncbi_table(self, fetched):
        dl.ensure_reference_data(has_ncbi_accessions=True)
        assert fetched == ["ncbi"]

    def test_an_ncbi_ref_tax_pull_needs_only_the_ncbi_table(self, fetched):
        dl.ensure_reference_data(wanted_ref_tax="Bacteria", source="ncbi")
        assert fetched == ["ncbi"]

    def test_a_gtdb_ref_tax_pull_needs_both_tables(self, fetched):
        """GTDB rows are resolved against NCBI assembly info, so both are required."""
        dl.ensure_reference_data(wanted_ref_tax="Bacteria", source="gtdb")
        assert fetched == ["ncbi", "gtdb"]

    @pytest.mark.parametrize("source", ["GTDB", "  gtdb  ", "Gtdb"])
    def test_the_source_is_matched_case_and_space_insensitively(self, fetched, source):
        dl.ensure_reference_data(wanted_ref_tax="Bacteria", source=source)
        assert fetched == ["ncbi", "gtdb"]

    def test_a_source_without_a_taxon_is_ignored(self, fetched):
        # source alone doesn't mean a taxonomy pull is happening
        dl.ensure_reference_data(source="gtdb")
        assert fetched == []
