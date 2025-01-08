import traitlets as tl

import numpy as np
from aiidalab_qe.common.mixins import HasInputStructure
from aiidalab_qe.common.panel import ConfigurationSettingsModel

from aiida_phonopy.data.preprocess import PreProcessData
from aiida.plugins import DataFactory
import sys
import os

HubbardStructureData = DataFactory("quantumespresso.hubbard_structure")
from aiida_vibroscopy.calculations.spectra_utils import get_supercells_for_hubbard
from aiida_vibroscopy.workflows.phonons.base import get_supercell_hubbard_structure

# spinner for waiting time (supercell estimations)
spinner_html = """
<style>
@keyframes spin {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(360deg); }
}

.spinner {
  display: inline-block;
  width: 15px;
  height: 15px;
}

.spinner div {
  width: 100%;
  height: 100%;
  border: 4px solid #f3f3f3;
  border-top: 4px solid #3498db;
  border-radius: 50%;
  animation: spin 1s linear infinite;
}
</style>
<div class="spinner">
  <div></div>
</div>
"""


def disable_print(func):
    def wrapper(*args, **kwargs):
        # Save the current standard output
        original_stdout = sys.stdout
        # Redirect standard output to os.devnull
        sys.stdout = open(os.devnull, "w")
        try:
            # Call the function
            result = func(*args, **kwargs)
        finally:
            # Restore the original standard output
            sys.stdout.close()
            sys.stdout = original_stdout
        return result

    return wrapper


class VibroConfigurationSettingsModel(ConfigurationSettingsModel, HasInputStructure):
    title = "Vibrational Settings"

    dependencies = [
        "input_structure",
    ]

    simulation_type_options = tl.List(
        trait=tl.List(tl.Union([tl.Unicode(), tl.Int()])),
        default_value=[
            ("IR/Raman, Phonon, Dielectric, INS properties", 1),
            ("IR/Raman and Dielectric in Primitive Cell Approach", 2),
            ("Phonons for non-polar materials and INS", 3),
            ("Dielectric properties", 4),
        ],
    )
    simulation_type = tl.Int(1)

    symmetry_symprec = tl.Float(1e-5)
    supercell_x = tl.Int(2)
    supercell_y = tl.Int(2)
    supercell_z = tl.Int(2)

    # Control for disable the supercell widget

    disable_x = tl.Bool(False)
    disable_y = tl.Bool(False)
    disable_z = tl.Bool(False)

    supercell = tl.List(
        trait=tl.Int(),
        default_value=[2, 2, 2],
    )
    supercell_number_estimator = tl.Unicode(
        "Click the button to estimate the supercell size."
    )

    def get_model_state(self):
        return {
            "simulation_type": self.simulation_type,
            "symmetry_symprec": self.symmetry_symprec,
            "supercell": self.supercell,
        }

    def set_model_state(self, parameters: dict):
        self.simulation_type = parameters.get("simulation_type", 1)
        self.symmetry_symprec = parameters.get("symmetry_symprec", 1e-5)
        self.supercell = parameters.get("supercell", [2, 2, 2])
        self.supercell_x, self.supercell_y, self.supercell_z = self.supercell

    def reset(self):
        with self.hold_trait_notifications():
            self.simulation_type = 1
            self.symmetry_symprec = self._get_default("symmetry_symprec")
            self.supercell = [2, 2, 2]
            self.supercell_x, self.supercell_y, self.supercell_z = self.supercell
            self.supercell_number_estimator = self._get_default(
                "supercell_number_estimator"
            )

    def _get_default(self, trait):
        return self._defaults.get(trait, self.traits()[trait].default_value)

    def on_input_structure_change(self, _=None):
        if not self.input_structure:
            self.reset()

        else:
            self.disable_x, self.disable_y, self.disable_z = True, True, True
            pbc = self.input_structure.pbc

            if pbc == (False, False, False):
                # No periodicity; fully disable and reset supercell
                self.supercell_x = self.supercell_y = self.supercell_z = 1
            elif pbc == (True, False, False):
                self.supercell_y = self.supercell_z = 1
                self.disable_x = False
                self.symmetry_symprec = 1e-3
            elif pbc == (True, True, False):
                self.supercell_z = 1
                self.disable_x = self.disable_y = False
            elif pbc == (True, True, True):
                self.disable_x = self.disable_y = self.disable_z = False

            self.supercell = [self.supercell_x, self.supercell_y, self.supercell_z]

    def suggest_supercell(self, _=None):
        """
        minimal supercell size for phonons, imposing a minimum lattice parameter of 15 A.
        """
        if self.input_structure and self.input_structure.pbc != (False, False, False):
            ase_structure = self.input_structure.get_ase()
            suggested_3D = 15 // np.array(ase_structure.cell.cellpar()[:3]) + 1

            # Update only dimensions that are not disabled
            if not self.disable_x:
                self.supercell_x = int(suggested_3D[0])
            if not self.disable_y:
                self.supercell_y = int(suggested_3D[1])
            if not self.disable_z:
                self.supercell_z = int(suggested_3D[2])

            # Sync the updated values to the supercell list
            self.supercell = [self.supercell_x, self.supercell_y, self.supercell_z]

        else:
            return

    def supercell_reset(self, _=None):
        if not self.disable_x:
            self.supercell_x = self._get_default("supercell_x")
        if not self.disable_y:
            self.supercell_y = self._get_default("supercell_x")
        if not self.disable_z:
            self.supercell_z = self._get_default("supercell_x")
        self.supercell = [self.supercell_x, self.supercell_y, self.supercell_z]

    def reset_symprec(self, _=None):
        self.symmetry_symprec = (
            self._get_default("symmetry_symprec")
            if self.input_structure.pbc != (True, False, False)
            else 1e-3
        )
        self.supercell_number_estimator = self._get_default(
            "supercell_number_estimator"
        )

    @disable_print
    def _estimate_supercells(self, _=None):
        if self.input_structure:
            self.supercell_number_estimator = spinner_html

            preprocess_data = PreProcessData(
                structure=self.input_structure,
                supercell_matrix=[
                    [self.supercell_x, 0, 0],
                    [0, self.supercell_y, 0],
                    [0, 0, self.supercell_z],
                ],
                symprec=self.symmetry_symprec,
                distinguish_kinds=False,
                is_symmetry=True,
            )

            if isinstance(self.input_structure, HubbardStructureData):
                supercell = get_supercell_hubbard_structure(
                    self.input_structure,
                    self.input_structure,
                    metadata={"store_provenance": False},
                )
                supercells = get_supercells_for_hubbard(
                    preprocess_data=preprocess_data,
                    ref_structure=supercell,
                    metadata={"store_provenance": False},
                )
            else:
                supercells = preprocess_data.get_supercells_with_displacements()

            self.supercell_number_estimator = f"{len(supercells)}"

        return
