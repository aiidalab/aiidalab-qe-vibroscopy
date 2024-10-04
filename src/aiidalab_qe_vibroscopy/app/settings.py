# -*- coding: utf-8 -*-
"""Panel for PhononWorkchain plugin.

Authors:

    * Miki Bonacci <miki.bonacci@psi.ch>
    Inspired by Xing Wang <xing.wang@psi.ch>
"""

import ipywidgets as ipw
import traitlets as tl
import numpy as np

from aiida import orm
from aiidalab_qe.common.panel import Panel


import sys
import os

from aiida.plugins import DataFactory

HubbardStructureData = DataFactory("quantumespresso.hubbard_structure")

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


class Setting(Panel):
    title = "Vibrational Settings"

    simulation_mode = [
        ("IR/Raman, Phonon, Dielectric, INS properties", 1),
        ("IR/Raman and Dielectric in Primitive Cell Approach", 2),
        ("Phonons for non-polar materials and INS", 3),
        ("Dielectric properties", 4),
    ]

    input_structure = tl.Instance(orm.StructureData, allow_none=True)

    def __init__(self, **kwargs):
        self.settings_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h4>Vibrational Settings</h4></div>"""
        )
        self.settings_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Calculations are performed using the <b><a href="https://aiida-vibroscopy.readthedocs.io/en/latest/"
        target="_blank">aiida-vibroscopy</b></a> plugin (L. Bastonero and N. Marzari, <a href="https://www.nature.com/articles/s41524-024-01236-3"
        target="_blank">npj Comput. Mater. <b>10</b>, 55, 2024</a>).
            The plugin employes the finite-displacement and finite-field approach. Raman spectra are simulated in the first-order non-resonant regime.
            </div>""",
            layout=ipw.Layout(width="400"),
        )

        self.use_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h5>Available vibrational properties:</h5></div>"""
        )

        self.use_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            <li style="margin-right: 10px; list-style-type: none; display: inline-block;">&#8226; <em>IR/Raman spectra</em>: both single crystal and powder samples.</li> <br>
            <li style="margin-right: 10px; list-style-type: none; display: inline-block;">&#8226; <em>Phonons properties</em>: bands, density of states and thermal properties (Helmoltz free energy, entropy and specific heat at constant volume).</li> <br>
            <li style="list-style-type: none; display: inline-block;">&#8226; <em>Dielectric properties</em>: Born charges, high-frequency dielectric tensor, non-linear optical susceptibility and raman tensors .</li> <br>
            <li style="list-style-type: none; display: inline-block;">&#8226; <em>Inelastic neutron scattering (INS)</em>: dynamic structure factor and powder intensity maps.</li> <br> <br>
            </div>""",
            layout=ipw.Layout(width="400"),
        )

        self.hint_button_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Select a supercell size for Phonon properties. Larger supercells increase computational costs. A 2x2x2 supercell is usually adequate.<br>
            Utilize the <em>Size hint</em> button for an estimate, maintaining a minimum lattice vector magnitude of 15Å along periodic directions.
            </div>""",
        )

        self.calc_options_description = ipw.HTML("Select calculation:")
        self.calc_options = ipw.Dropdown(
            options=self.simulation_mode,
            layout=ipw.Layout(width="450px"),
            value=self.simulation_mode[0][1],
        )

        self.calc_options.observe(self._display_supercell, names="value")

        # start Supercell

        self.supercell = [2, 2, 2]

        def change_supercell(_=None):
            self.supercell = [
                self._sc_x.value,
                self._sc_y.value,
                self._sc_z.value,
            ]

        if self.input_structure:
            pbc = self.input_structure.pbc
        else:
            pbc = (True, True, True)

        for elem, periodic in zip(["x", "y", "z"], pbc):
            # periodic allows support of hints also for 2D, 1D.
            setattr(
                self,
                "_sc_" + elem,
                ipw.BoundedIntText(
                    value=2 if periodic else 1,
                    min=1,
                    layout={"width": "40px"},
                    disabled=False if periodic else True,
                ),
            )
        for elem in [self._sc_x, self._sc_y, self._sc_z]:
            elem.observe(change_supercell, names="value")
            elem.observe(self._activate_estimate_supercells, names="value")

        self.supercell_selector = ipw.HBox(
            children=[
                ipw.HTML(
                    description="Supercell size:",
                    style={"description_width": "initial"},
                )
            ]
            + [
                self._sc_x,
                self._sc_y,
                self._sc_z,
            ],
        )

        ## start supercell hint:

        # supercell data
        self.supercell_hint_button = ipw.Button(
            description="Size hint",
            disabled=False,
            layout=ipw.Layout(width="100px"),
            button_style="info",
        )
        # supercell hint (15A lattice params)
        self.supercell_hint_button.on_click(self._suggest_supercell)

        # reset supercell
        self.supercell_reset_button = ipw.Button(
            description="Reset hint",
            disabled=False,
            layout=ipw.Layout(width="100px"),
            button_style="warning",
        )
        # supercell reset reaction
        self.supercell_reset_button.on_click(self._reset_supercell)

        # Estimate supercell button
        self.supercell_estimate_button = ipw.Button(
            description="Estimate number of supercells ➡",
            disabled=False,
            layout=ipw.Layout(width="240px", display="none"),
            button_style="info",
            tooltip="Number of supercells for phonons calculations;\nwarning: for large systems, this may take some time.",
        )
        # supercell reset reaction
        self.supercell_estimate_button.on_click(self._estimate_supercells)

        # Estimate the number of supercells for frozen phonons.
        self.supercell_number_estimator = ipw.HTML(
            # description="Number of supercells:",
            value="?",
            style={"description_width": "initial"},
            layout=ipw.Layout(display="none"),
        )

        ## end supercell hint.

        self.supercell_widget = ipw.VBox(
            [
                self.hint_button_help,
                ipw.HBox(
                    [
                        self.supercell_selector,
                        self.supercell_hint_button,
                        self.supercell_reset_button,
                        self.supercell_estimate_button,  # I do it on request, as it can take long time.
                        self.supercell_number_estimator,
                    ],
                ),
            ]
        )
        self.supercell_widget.layout.display = "block"
        # end Supercell.

        self.symmetry_symprec = ipw.FloatText(
            value=1e-5,
            max=1,
            min=1e-7,  # Ensure the value is always positive
            step=1e-4,  # Step value of 1e-4
            description="Symmetry tolerance (symprec):",
            style={"description_width": "initial"},
            layout={"width": "300px"},
        )
        self.symmetry_symprec.observe(self._activate_estimate_supercells, "value")

        # reset supercell
        self.symmetry_symprec_reset_button = ipw.Button(
            description="Reset symprec",
            disabled=False,
            layout=ipw.Layout(width="125px"),
            button_style="warning",
        )
        # supercell reset reaction
        self.symmetry_symprec_reset_button.on_click(self._reset_symprec)

        self.children = [
            ipw.VBox(
                [
                    ipw.VBox(
                        [
                            self.settings_title,
                            self.settings_help,
                        ]
                    ),
                    ipw.VBox(
                        [
                            self.use_title,
                            self.use_help,
                        ]
                    ),
                ]
            ),
            ipw.HBox(
                [
                    self.calc_options_description,
                    self.calc_options,
                ],
            ),
            self.supercell_widget,
            ipw.HBox(
                [
                    self.symmetry_symprec,
                    self.symmetry_symprec_reset_button,
                ],
            ),
        ]

        super().__init__(**kwargs)

        # we define a block for the estimation of the supercell if we ask for hint,
        # so that we call the estimator only at the end of the supercell hint generator,
        # and now each time after the x, y, z generation (i.e., we don't lose time).
        # see the methods below.
        self.block = False

    @tl.observe("input_structure")
    def _update_input_structure(self, change):
        if self.input_structure:
            for direction, periodic in zip(
                [self._sc_x, self._sc_y, self._sc_z], self.input_structure.pbc
            ):
                direction.value = 2 if periodic else 1
                direction.disabled = False if periodic else True

            self.supercell_number_estimator.layout.display = (
                "block" if len(self.input_structure.sites) <= 30 else "none"
            )
            self.supercell_estimate_button.layout.display = (
                "block" if len(self.input_structure.sites) <= 30 else "none"
            )
        else:
            self.supercell_number_estimator.layout.display = "none"
            self.supercell_estimate_button.layout.display = "none"

    def _display_supercell(self, change):
        selected = change["new"]
        if selected in [1, 3]:
            self.supercell_widget.layout.display = "block"
        else:
            self.supercell_widget.layout.display = "none"

    def _suggest_supercell(self, _=None):
        """
        minimal supercell size for phonons, imposing a minimum lattice parameter of 15 A.
        """
        if self.input_structure:
            s = self.input_structure.get_ase()
            suggested_3D = 15 // np.array(s.cell.cellpar()[:3]) + 1

            # if disabled, it means that it is a non-periodic direction.
            # here we manually unobserve the `_activate_estimate_supercells`, so it is faster
            # and only compute when all the three directions are updated
            self.block = True
            for direction, suggested, original in zip(
                [self._sc_x, self._sc_y, self._sc_z], suggested_3D, s.cell.cellpar()[:3]
            ):
                direction.value = suggested if not direction.disabled else 1
            self.block = False
            self._activate_estimate_supercells()
        else:
            return

    def _activate_estimate_supercells(self, _=None):
        self.supercell_estimate_button.disabled = False
        self.supercell_number_estimator.value = "?"

    # @tl.observe("input_structure")
    @disable_print
    def _estimate_supercells(self, _=None):
        """_summary_

        Estimate the number of supercells to be computed for frozen phonon calculation.
        """
        if self.block:
            return

        symprec_value = self.symmetry_symprec.value

        self.symmetry_symprec.value = max(1e-5, min(symprec_value, 1))

        self.supercell_number_estimator.value = spinner_html

        from aiida_phonopy.data.preprocess import PreProcessData

        if self.input_structure:
            preprocess_data = PreProcessData(
                structure=self.input_structure,
                supercell_matrix=[
                    [self._sc_x.value, 0, 0],
                    [0, self._sc_y.value, 0],
                    [0, 0, self._sc_z.value],
                ],
                symprec=self.symmetry_symprec.value,
                distinguish_kinds=False,
                is_symmetry=True,
            )

            supercells = preprocess_data.get_supercells_with_displacements()

            # for now, we comment the following part, as the HubbardSD is generated in the submission step.
            """if isinstance(self.input_structure, HubbardStructureData):
                from aiida_vibroscopy.calculations.spectra_utils import get_supercells_for_hubbard
                from aiida_vibroscopy.workflows.phonons.base import get_supercell_hubbard_structure
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
            """
            self.supercell_number_estimator.value = f"{len(supercells)}"
            self.supercell_estimate_button.disabled = True

        return

    def _reset_supercell(self, _=None):
        if self.input_structure is not None:
            reset_supercell = []
            self.block = True
            for direction, periodic in zip(
                [self._sc_x, self._sc_y, self._sc_z], self.input_structure.pbc
            ):
                reset_supercell.append(2 if periodic else 1)
            (self._sc_x.value, self._sc_y.value, self._sc_z.value) = reset_supercell
            self.block = False
            self._activate_estimate_supercells()
        return

    def _reset_symprec(self, _=None):
        self.symmetry_symprec.value = 1e-5
        self._activate_estimate_supercells()
        return

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        return {
            "simulation_mode": self.calc_options.value,
            "supercell_selector": self.supercell,
            "symmetry_symprec": self.symmetry_symprec.value,
        }

    def set_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.calc_options.value = input_dict.get("simulation_mode", 1)
        self.supercell = input_dict.get("supercell_selector", [2, 2, 2])
        self.symmetry_symprec.value = input_dict.get("symmetry_symprec", 1e-5)
        self._sc_x.value, self._sc_y.value, self._sc_z.value = self.supercell

    def reset(self):
        """Reset the panel"""
        self.calc_options.value = 1
        self.supercell = [2, 2, 2]
        self.symmetry_symprec.value = 1e-5
        self._sc_x.value, self._sc_y.value, self._sc_z.value = self.supercell
