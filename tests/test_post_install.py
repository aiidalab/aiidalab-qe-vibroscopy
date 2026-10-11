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
    monkeypatch.setattr(plugin_cli, "machine", "x86_64")

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


def test_verify_post_install_checks_arm_dependency_versions(monkeypatch):
    monkeypatch.setattr(plugin_cli, "load_profile", lambda _profile: None)
    monkeypatch.setattr(plugin_cli, "load_code", lambda _label: object())
    monkeypatch.setattr(plugin_cli, "machine", "aarch64")
    monkeypatch.setattr(
        plugin_cli.metadata,
        "version",
        lambda package: {"scipy": "1.13.1", "spglib": "2.5.0"}[package],
    )

    result = CliRunner().invoke(plugin_cli.cli, ["verify-post-install"])

    assert result.exit_code == 0
