from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from aiidalab_qe_vibroscopy import __main__ as plugin_cli


def test_verify_post_install_passes(monkeypatch, tmp_path):
    executable = tmp_path / "phonopy"
    executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setattr(plugin_cli, "load_profile", lambda _profile: None)
    monkeypatch.setattr(
        plugin_cli,
        "load_code",
        lambda _label: SimpleNamespace(filepath_executable=str(executable)),
    )

    result = CliRunner().invoke(plugin_cli.cli, ["verify-post-install"])

    assert result.exit_code == 0
    assert "setup is valid" in result.output


def test_verify_post_install_fails_when_code_is_missing(monkeypatch):
    from aiida.common.exceptions import NotExistent

    monkeypatch.setattr(plugin_cli, "load_profile", lambda _profile: None)

    def missing_code(_label):
        raise NotExistent("missing")

    monkeypatch.setattr(plugin_cli, "load_code", missing_code)

    result = CliRunner().invoke(plugin_cli.cli, ["verify-post-install"])

    assert result.exit_code != 0
    assert "phonopy@localhost is not configured" in result.output


def test_setup_phonopy_only_registers_code(monkeypatch):
    from aiida.common.exceptions import NotExistent

    commands = []
    monkeypatch.setattr(plugin_cli, "load_profile", lambda: None)
    monkeypatch.setattr(
        plugin_cli,
        "load_code",
        lambda _label: (_ for _ in ()).throw(NotExistent("missing")),
    )
    monkeypatch.setattr(plugin_cli.shutil, "which", lambda _name: "/usr/bin/phonopy")
    monkeypatch.setattr(
        plugin_cli.subprocess,
        "run",
        lambda command, **_kwargs: (
            commands.append(command)
            or plugin_cli.subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    result = CliRunner().invoke(plugin_cli.cli, ["setup-phonopy"])

    assert result.exit_code == 0, result.output
    assert commands == [
        [
            "verdi",
            "code",
            "create",
            "core.code.installed",
            "--non-interactive",
            "--label",
            "phonopy",
            "--default-calc-job-plugin",
            "phonopy.phonopy",
            "--computer",
            "localhost",
            "--filepath-executable",
            "/usr/bin/phonopy",
        ]
    ]


@pytest.mark.parametrize(
    ("path_state", "expected_message"),
    [
        ("unset", "phonopy@localhost executable is not set"),
        ("missing", "phonopy@localhost executable is missing"),
        ("directory", "phonopy@localhost executable is not a file"),
        ("not_executable", "phonopy@localhost executable is not executable"),
    ],
)
def test_verify_post_install_reports_executable_path_state(
    monkeypatch, tmp_path, path_state, expected_message
):
    executable = tmp_path / "phonopy"
    if path_state == "unset":
        filepath = None
    elif path_state == "directory":
        executable.mkdir()
        filepath = str(executable)
    elif path_state == "not_executable":
        executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        executable.chmod(0o644)
        filepath = str(executable)
    else:
        filepath = str(executable)
    monkeypatch.setattr(plugin_cli, "load_profile", lambda _profile: None)
    monkeypatch.setattr(
        plugin_cli,
        "load_code",
        lambda _label: SimpleNamespace(filepath_executable=filepath),
    )

    result = CliRunner().invoke(plugin_cli.cli, ["verify-post-install"])

    assert result.exit_code != 0
    assert expected_message in result.output


def test_verify_post_install_runs_registered_executable(monkeypatch, tmp_path):
    executable = tmp_path / "phonopy"
    executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setattr(plugin_cli, "load_profile", lambda _profile: None)
    monkeypatch.setattr(
        plugin_cli,
        "load_code",
        lambda _label: SimpleNamespace(filepath_executable=str(executable)),
    )
    commands = []

    def run(command, **kwargs):
        commands.append(command)
        assert kwargs["timeout"] == 30
        return plugin_cli.subprocess.CompletedProcess(command, 0, "Phonopy help", "")

    monkeypatch.setattr(plugin_cli.subprocess, "run", run)

    result = CliRunner().invoke(plugin_cli.cli, ["verify-post-install"])

    assert result.exit_code == 0
    assert commands == [[str(executable), "--help"]]


@pytest.mark.parametrize(
    "error",
    [OSError("cannot execute"), plugin_cli.subprocess.TimeoutExpired("phonopy", 30)],
)
def test_verify_post_install_reports_startup_error(monkeypatch, tmp_path, error):
    executable = tmp_path / "phonopy"
    executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setattr(plugin_cli, "load_profile", lambda _profile: None)
    monkeypatch.setattr(
        plugin_cli,
        "load_code",
        lambda _label: SimpleNamespace(filepath_executable=str(executable)),
    )

    def run(*_args, **_kwargs):
        raise error

    monkeypatch.setattr(plugin_cli.subprocess, "run", run)

    result = CliRunner().invoke(plugin_cli.cli, ["verify-post-install"])

    assert result.exit_code != 0
    assert "Could not start AiiDA code phonopy@localhost" in result.output


def test_verify_post_install_reports_broken_executable(monkeypatch, tmp_path):
    executable = tmp_path / "phonopy"
    executable.write_text(
        "#!/bin/sh\necho 'ImportError: broken native library' >&2\nexit 1\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    monkeypatch.setattr(plugin_cli, "load_profile", lambda _profile: None)
    monkeypatch.setattr(
        plugin_cli,
        "load_code",
        lambda _label: SimpleNamespace(filepath_executable=str(executable)),
    )

    result = CliRunner().invoke(plugin_cli.cli, ["verify-post-install"])

    assert result.exit_code != 0
    assert "exit code 1" in result.output
    assert "ImportError: broken native library" in result.output


def test_setup_phonopy_preserves_existing_code(monkeypatch):
    monkeypatch.setattr(plugin_cli, "load_profile", lambda: None)
    monkeypatch.setattr(plugin_cli, "load_code", lambda _label: SimpleNamespace())
    commands = []
    monkeypatch.setattr(
        plugin_cli.subprocess,
        "run",
        lambda command, **_kwargs: commands.append(command),
    )

    result = CliRunner().invoke(plugin_cli.cli, ["setup-phonopy"])

    assert result.exit_code == 0
    assert "already installed" in result.output
    assert commands == []
