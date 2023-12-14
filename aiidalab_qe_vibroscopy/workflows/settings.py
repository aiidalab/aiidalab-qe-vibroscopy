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

    def __init__(self, **kwargs):
        self.settings_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h4>Settings</h4></div>"""
        )
        self.settings_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Calculations are performed using the <b><a href="https://aiida-vibroscopy.readthedocs.io/en/latest/"
        target="_blank">aiida-vibroscopy</b></a> plugin (L. Bastonero and N. Marzari, <a href="https://arxiv.org/abs/2308.04308"
        target="_blank">Automated all-functionals infrared and Raman spectra</a>). <br>
            The plugin employes the finite-displacement and finite-field approach, for Phonon dispersion please select a supercell size:
            the larger the supercell, the larger the computational cost of the simulations. Usually, a 2x2x2 supercell should be enough. <br>
            </div>"""
        )

        self.use_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            The plugin is capable to run the following calculations: <br>
            <li style="margin-right: 10px; list-style-type: none; display: inline-block;">&#8226; IR/Raman spectrum.</li>
            <li style="margin-right: 10px; list-style-type: none; display: inline-block;">&#8226; Phonon band dispersion.</li>
            <li style="list-style-type: none; display: inline-block;">&#8226; Phonon Projected Density of States.</li>
            <li style="list-style-type: none; display: inline-block;">&#8226; Thermal Properties.</li>
            <li style="list-style-type: none; display: inline-block;">&#8226; Dielectric Properties.</li>
            </div>"""
        )

        self.polar_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 5px; padding-bottom: 5px">
            If the material is polar, more
            accurate phonon properties interpolation is performed.
            </div>"""
        )

        self.workchain_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
            value="moderate",
        )

        self.calc_options = ipw.Dropdown(
            description="Calculation:",
            options=[
                ["Raman spectrum", "raman"],
                ["IR spectrum", "ir"],
                ["Phonons Bands", "ph_bands"],
                ["Phonons Pdos", "ph_pdos"],
                ["Thermal properties", "ph_therm"],
                ["Dielectric properties", "dielectric"],
            ],
            value="raman",
            style={"description_width": "initial"},
        )

        self.supercell_out = ipw.Output()

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

        # to trigger Dielectric property = Raman... FOR POLAR MATERIALS.
        self.material_is_polar_ = ipw.ToggleButtons(
            options=[("Off", "off"), ("On", "on")],
            value="off",
            style={"description_width": "initial"},
        )
        # self.material_is_polar_.observe(self._onclick_material_is_polar, "value")
        self.calc_options.observe(self._display_supercell, names="value")
        self.children = [
            self.settings_title,
            self.settings_help,
            self.use_help,
            ipw.HBox(
                [
                    self.calc_options,
                    self.supercell_out,
                ]
            ),
            self.polar_help,
            ipw.HBox(
                children=[
                    ipw.Label(
                        "Material is polar:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.material_is_polar_,
                ]
            ),
        ]
        super().__init__(**kwargs)

    def _display_supercell(self, change):
        selected = change["new"]
        if selected in ["ph_bands", "ph_pdos", "ph_therm"]:
            with self.supercell_out:
                clear_output()
                display(self.supercell_widget)
        else:
            with self.supercell_out:
                clear_output()

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        if isinstance(self.calc_options, str):
            return {
                "calc_options": self.calc_options,
                "material_is_polar": self.material_is_polar_.value,
                "supercell_selector": self.supercell,
            }
        return {
            "calc_options": self.calc_options.value,
            "material_is_polar": self.material_is_polar_.value,
            "supercell_selector": self.supercell,
        }

    def load_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.material_is_polar_.value = input_dict.get("material_is_polar", "off")
        self.calc_options.value = input_dict.get("calc_options", "raman")
        self.supercell = input_dict.get("supercell_selector", [2, 2, 2])

    def reset(self):
        """Reset the panel"""
        self.material_is_polar_.value = "off"
        if isinstance(self.calc_options, str):
            self.calc_options = "raman"
        else:
            self.calc_options.value = "raman"
        self.supercell = [2, 2, 2]
