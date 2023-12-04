from aiida.common.exceptions import NotExistent
from subprocess import run
from aiida.orm import load_code
from aiida import load_profile


def install_phonopy():
    load_profile()
    try:
        load_code("phonopy@localhost")
    except NotExistent:
        run(
            [
                "verdi",
                "code",
                "create",
                "core.code.installed",
                "--non-interactive",
                "--label",
                "phonopy",
                "--description",
                "phonopy setup by AiiDAlab.",
                "--default-calc-job-plugin",
                "phonopy.phonopy",
                "--computer",
                "localhost",
                "--filepath-executable",
                "/opt/conda/bin/phonopy",
            ],
            check=True,
            capture_output=True,
        )
