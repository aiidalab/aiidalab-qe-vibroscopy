from aiida.common.exceptions import NotExistent
import subprocess
from aiida.orm import load_code
from aiida import load_profile
import shutil

"""
Automatic installation of the phonopy code.

(1) if not already installed: pip install phonopy --user
(2) if it doesn't already exist, create a symbolic link for Phonopy: ln -s <phonopy path> /opt/conda/bin/phonopy (try which phonopy to see the phonopy path)
(3) invoke phonopy_setup
"""


def install_phonopy():
    load_profile()
    try:
        load_code("phonopy@localhost")
    except NotExistent:
        # Use shutil.which to find the path of the phonopy executable
        phonopy_path = shutil.which("phonopy")
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
        raise Warning("Code phonopy@localhost already installed!")


# Called when the script is run directly
if __name__ == "__main__":
    install_phonopy()
