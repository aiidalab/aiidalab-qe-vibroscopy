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
from IPython.display import clear_output, display


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
            width="500px",
            button_style="info",
        )
        # supercell hint (15A lattice params)
        self.supercell_hint_button.on_click(self._suggest_supercell)

        ## end supercell hint.

        self.supercell_widget = ipw.VBox(
            [
                self.hint_button_help,
                ipw.HBox(
                    [self.supercell_selector, self.supercell_hint_button],
                ),
            ]
        )
        self.supercell_widget.layout.display = "block"
        # end Supercell.

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
        ]

        super().__init__(**kwargs)

    @tl.observe("input_structure")
    def _update_input_structure(self, change):
        if self.input_structure is not None:
            for direction, periodic in zip(
                [self._sc_x, self._sc_y, self._sc_z], self.input_structure.pbc
            ):
                direction.value = 2 if periodic else 1
                direction.disabled = False if periodic else True

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
            for direction, suggested, original in zip(
                [self._sc_x, self._sc_y, self._sc_z], suggested_3D, s.cell.cellpar()[:3]
            ):
                direction.value = suggested if not direction.disabled else 1
        else:
            return

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        return {
            "simulation_mode": self.calc_options.value,
            "supercell_selector": self.supercell,
        }

    def load_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.calc_options.value = input_dict.get("simulation_mode", 1)
        self.supercell = input_dict.get("supercell_selector", [2, 2, 2])

    def reset(self):
        """Reset the panel"""
        self.calc_options.value = 1
        self.supercell = [2, 2, 2]
