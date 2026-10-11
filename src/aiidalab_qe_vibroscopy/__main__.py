import os
import shutil
import subprocess
from pathlib import Path

import click
from aiida import load_profile
from aiida.common.exceptions import NotExistent
from aiida.orm import load_code

"""
Automatic installation of the phonopy code.

we suppose phonopy is already installed (pip install phonopy --user).
So we only setup in AiiDA.
"""

import platform

try:
    machine = platform.machine()
except Exception:
    machine = "aarch"


@click.group()
def cli():
    pass


@cli.command(help="Setup phonopy@localhost in the current AiiDA database.")
def setup_phonopy():
    load_profile()
    try:
        load_code("phonopy@localhost")
    except NotExistent:
        # Use shutil.which to find the path of the phonopy executable
        phonopy_path = shutil.which("phonopy")
        if not phonopy_path:
            raise FileNotFoundError(
                "Phonopy code is not found in PATH.  \
                                    You should update your PATH. If you have not phonopy in , \
                                    your environment, install the code via \
                                    `pip install phonopy --user`."
            )
        # Construct the command as a list of arguments
        command = [
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
            phonopy_path,
        ]

        # Use subprocess.run to run the command
        subprocess.run(command, check=True)
    else:
        print("Code phonopy@localhost is already installed! Nothing to do here.")

    if (
        "arm" in machine
    ):  # we need to install scipy from conda only if arm machine. OPENBLAS broken in the pip one.
        command = [
            "conda",
            "install",
            "scipy==1.13.1",
            "--y",
        ]
        # Use subprocess.run to run the command
        subprocess.run(command, check=True)

        command = [
            "pip",
            "install",
            "spglib==2.5.0",
            "--y",
        ]
        # Use subprocess.run to run the command
        subprocess.run(command, check=True)


@cli.command(help="Verify the phonopy post-install setup.")
@click.option("-p", "--profile", default="default", help="AiiDA profile name.")
def verify_post_install(profile):
    from importlib_metadata import PackageNotFoundError, version

    load_profile(profile)
    try:
        code = load_code("phonopy@localhost")
    except NotExistent as error:
        raise click.ClickException(
            "AiiDA code phonopy@localhost is not configured. Run setup-phonopy."
        ) from error

    executable = code.filepath_executable
    if not executable:
        raise click.ClickException("AiiDA code phonopy@localhost executable is not set")

    filepath = Path(executable)
    failure_prefix = "AiiDA code phonopy@localhost executable is"
    if not filepath.exists():
        raise click.ClickException(f"{failure_prefix} missing: {executable}")
    if not filepath.is_file():
        raise click.ClickException(f"{failure_prefix} not a file: {executable}")
    if not os.access(filepath, os.X_OK):
        raise click.ClickException(f"{failure_prefix} not executable: {executable}")

    if "arm" in machine:
        expected_versions = {"scipy": "1.13.1", "spglib": "2.5.0"}
        for package, expected_version in expected_versions.items():
            try:
                installed_version = version(package)
            except PackageNotFoundError as error:
                raise click.ClickException(
                    f"{package} is not installed; expected {expected_version}."
                ) from error
            if installed_version != expected_version:
                raise click.ClickException(
                    f"{package} {installed_version} is installed; "
                    f"expected {expected_version}."
                )

    click.echo("Vibroscopy post-install setup is valid.")


if __name__ == "__main__":
    cli()
