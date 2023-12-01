import pathlib
from setuptools import setup, find_packages
from aiida.common.exceptions import NotExistent
from subprocess import run
from aiida.orm import load_code
from aiida import load_profile
from setuptools.command.install import install

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()


class CustomInstallCommand(install):
    def run(self):
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
        else:
            raise RuntimeError(f"Code phonopy is already set up!")
        super().run()


setup(
    name="aiidalab-qe-vibroscopy",
    version="0.0.1",
    description="A aiidalab-qe plugin to calculate the vibrational properties of materials via aiida-vibroscopy package.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/mikibonacci/aiidalab-qe-vibroscopy",
    author="Miki Bonacci",
    author_email="mikibonacci@hotmail.it",
    classifiers=[],
    packages=find_packages(),
    entry_points={
        "aiidalab_qe.properties": [
            # "harmonic = aiidalab_qe_vibroscopy.harmonic:property",
            # "iraman = aiidalab_qe_vibroscopy.raman:property",
            # "dielectric = aiidalab_qe_vibroscopy.dielectric:property",
            # "phonons = aiidalab_qe_vibroscopy.phonons:property",
            "vibronic = aiidalab_qe_vibroscopy.workflows:property",
        ],
        "aiida.workflows": [
            "vibroscopy_app.vibro = aiidalab_qe_vibroscopy.workflows.vibroworkchain:VibroWorkChain",
        ],
    },
    install_requires=[
        "aiida-vibroscopy>=1.0.2",
        "aiida-phonopy>=1.1.3",
        "phonopy",
        "pre-commit",
    ],
    cmdclass={
        "install": CustomInstallCommand,
    },
    package_data={},
    python_requires=">=3.6",
)
