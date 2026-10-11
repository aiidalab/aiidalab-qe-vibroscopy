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


@pytest.mark.parametrize("machine", ["arm64", "aarch64"])
def test_verify_post_install_checks_arm_dependency_versions(
    monkeypatch, tmp_path, machine
):
    executable = tmp_path / "phonopy"
    executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setattr(plugin_cli, "load_profile", lambda _profile: None)
    monkeypatch.setattr(
        plugin_cli,
        "load_code",
        lambda _label: SimpleNamespace(filepath_executable=str(executable)),
    )
    monkeypatch.setattr(plugin_cli, "machine", machine)
    checked_packages = []

    def get_version(package):
        checked_packages.append(package)
        return {"scipy": "1.13.1", "spglib": "2.5.0"}[package]

    monkeypatch.setattr(
        plugin_cli.importlib_metadata,
        "version",
        get_version,
    )

    result = CliRunner().invoke(plugin_cli.cli, ["verify-post-install"])

    assert result.exit_code == 0
    assert set(checked_packages) == {"scipy", "spglib"}


@pytest.mark.parametrize("machine", ["arm64", "aarch64"])
@pytest.mark.parametrize(
    ("state", "package", "expected_message"),
    [
        ("missing", "scipy", "scipy is not installed; expected 1.13.1"),
        ("wrong", "scipy", "scipy 1.13.0 is installed; expected 1.13.1"),
        ("missing", "spglib", "spglib is not installed; expected 2.5.0"),
        ("wrong", "spglib", "spglib 2.4.0 is installed; expected 2.5.0"),
    ],
)
def test_verify_post_install_reports_arm_dependency_problems(
    monkeypatch, tmp_path, machine, state, package, expected_message
):
    executable = tmp_path / "phonopy"
    executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setattr(plugin_cli, "load_profile", lambda _profile: None)
    monkeypatch.setattr(
        plugin_cli,
        "load_code",
        lambda _label: SimpleNamespace(filepath_executable=str(executable)),
    )
    monkeypatch.setattr(plugin_cli, "machine", machine)

    def get_version(name):
        if name == package and state == "missing":
            raise plugin_cli.importlib_metadata.PackageNotFoundError(name)
        if name == "scipy":
            return "1.13.0" if package == name and state == "wrong" else "1.13.1"
        return "2.4.0" if package == name and state == "wrong" else "2.5.0"

    monkeypatch.setattr(plugin_cli.importlib_metadata, "version", get_version)

    result = CliRunner().invoke(plugin_cli.cli, ["verify-post-install"])

    assert result.exit_code != 0
    assert expected_message in result.output
