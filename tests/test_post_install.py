from click.testing import CliRunner

from aiidalab_qe_vibroscopy import __main__ as plugin_cli


def test_verify_post_install_passes(monkeypatch):
    monkeypatch.setattr(plugin_cli, "load_profile", lambda _profile: None)
    monkeypatch.setattr(plugin_cli, "load_code", lambda _label: object())
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
