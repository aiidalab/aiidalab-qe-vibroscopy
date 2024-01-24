# -*- coding: utf-8 -*-
"""Panel for PhononWorkchain plugin.

Authors:

    * Miki Bonacci <miki.bonacci@psi.ch>
    Inspired by Xing Wang <xing.wang@psi.ch>
"""
import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.common.panel import Panel
from IPython.display import clear_output, display


class Setting(Panel):
    title = "Vibrational Settings"
    choices = [
        (
            "Full: IR/Raman spectra, Phonon properties, Dielectric properties and Inelastic neutron scattering",
            1,
        ),
        ("IR/Raman and Dielectric properties in primitive cell approach", 2),
        (
            "Phonons and Inelastic neutron scattering without non-analytic corrections for force constants (non-polar materials)",
            3,
        ),
        ("Dielectric properties only", 4),
    ]

    def __init__(self, **kwargs):
        self.settings_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h4>Vibrational Settings</h4></div>"""
        )
        self.settings_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Calculations are performed using the <b><a href="https://aiida-vibroscopy.readthedocs.io/en/latest/"
        target="_blank">aiida-vibroscopy</b></a> plugin (L. Bastonero and N. Marzari, <a href="https://arxiv.org/abs/2308.04308"
        target="_blank">Automated all-functionals infrared and Raman spectra</a>).
            The plugin employes the finite-displacement and finite-field approach.
            </div>""",
            layout=ipw.Layout(width="400"),
        )

        self.use_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h5>Available vibrational properties:</h5></div>"""
        )

        self.use_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            The plugin is capable of simulating the following properties: <br>
            <li style="margin-right: 10px; list-style-type: none; display: inline-block;">&#8226; IR/Raman spectra, both single crystal and powder samples.</li> <br>
            <li style="margin-right: 10px; list-style-type: none; display: inline-block;">&#8226; Phonons properties: bands, density of states and thermal properties.</li> <br>
            <li style="list-style-type: none; display: inline-block;">&#8226; Dielectric properties: Born charges, high-frequency dielectric tensor, non-linear optical susceptibility and raman tensors .</li> <br>
            <li style="list-style-type: none; display: inline-block;">&#8226; Inelastic neutron scattering: dynamic structure factor and powder intensity maps.</li> <br> <br>
            For Phonon properties, please select a supercell size: the larger the supercell, the larger the computational cost of the simulations. Usually, a 2x2x2 supercell should be enough. <br>
            Raman spectra are simulated in the first-order non-resonant regime.
            </div>""",
            layout=ipw.Layout(width="400"),
        )

        self.workchain_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
            value="moderate",
        )

        self.calc_options = ipw.Dropdown(
            description="Select calculation:",
            options=self.choices,
            layout=ipw.Layout(width="750px"),
            value=1,
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

        for elem in ["x", "y", "z"]:
            setattr(
                self,
                "_sc_" + elem,
                ipw.BoundedIntText(
                    value=2, min=1, layout={"width": "40px"}, disabled=False
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

        self.supercell_widget = ipw.HBox(
            [self.supercell_selector],
            layout=ipw.Layout(justify_content="flex-start"),
        )
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
            self.calc_options,
            self.supercell_widget,
        ]

        super().__init__(**kwargs)

    def _display_supercell(self, change):
        selected = change["new"]
        if selected in [1, 3]:
            self.supercell_widget.layout.display = "block"
        else:
            self.supercell_widget.layout.display = "none"

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        if isinstance(self.calc_options, str):
            return {
                "simulation_mode": self.calc_options.value,
                "supercell_selector": self.supercell,
            }
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
