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


@cli.command(help="Verify the phonopy post-install setup.")
@click.option("-p", "--profile", default="default", help="AiiDA profile name.")
def verify_post_install(profile):
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

    try:
        result = subprocess.run(
            [executable, "--help"],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise click.ClickException(
            f"Could not start AiiDA code phonopy@localhost: {error}"
        ) from error

    if result.returncode != 0:
        details = (result.stderr or result.stdout).strip()[:500]
        raise click.ClickException(
            f"AiiDA code phonopy@localhost failed to start "
            f"(exit code {result.returncode}): {details}"
        )

    click.echo("Vibroscopy post-install setup is valid.")


if __name__ == "__main__":
    cli()
