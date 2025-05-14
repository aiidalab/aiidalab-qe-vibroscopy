from aiida.common.exceptions import NotExistent
import subprocess
from aiida.orm import load_code
from aiida import load_profile
import shutil
import click

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


if __name__ == "__main__":
    cli()
