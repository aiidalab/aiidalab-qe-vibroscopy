# -*- coding: utf-8 -*-
"""Panel for PhononWorkchain plugin.

Authors:

    * Miki Bonacci <miki.bonacci@psi.ch>
    Inspired by Xing Wang <xing.wang@psi.ch>
"""
import ipywidgets as ipw
from aiida.orm import Float, Int, Str

from aiidalab_qe.common.panel import Panel


class Setting(Panel):
    title = "Phonons Settings"
    identifier = "harmonic"

    def __init__(self, **kwargs):
        self.settings_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h4>Phonons settings</h4></div>"""
        )
        self.settings_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Please select the phonon-related properties to be computed in this simulation.
            You can select also the size of the supercell to be used; usually, a 2X2X2 supercell is enough
            to obtain converged results.
            </div>"""
        )

        self.polar_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 5px; padding-bottom: 5px">
            If the material is polar, also 3rd order derivatives will be computed and more
            accurate phonon band interpolation is performed.
            </div>"""
        )

        self.workchain_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
            value="moderate",
        )

        # I want to be able to select more than only one... this has to change at the PhononWorkChain level.
        self.phonon_property = ipw.Dropdown(
            options=[
                ["band structure", "BANDS"],
                ["density of states (DOS)", "DOS"],
                ["thermal properties", "THERMODYNAMIC"],
                ["force constants", "NONE"],
                ["none", "none"],
            ],
            value="BANDS",
            description="Phonon property:",
            disabled=False,
            style={"description_width": "initial"},
        )

        # 1. Supercell
        self.supercell = [1, 1, 1]

        def change_supercell(_=None):
            self.supercell = [
                _supercell[0].value,
                _supercell[1].value,
                _supercell[2].value,
            ]

        _supercell = [
            ipw.BoundedIntText(value=1, min=1, layout={"width": "40px"}),
            ipw.BoundedIntText(value=1, min=1, layout={"width": "40px"}),
            ipw.BoundedIntText(value=1, min=1, layout={"width": "40px"}),
        ]
        for elem in _supercell:
            elem.observe(change_supercell, names="value")
        self.supercell_selector = ipw.HBox(
            children=[
                ipw.HTML(
                    description="Supercell size:",
                    style={"description_width": "initial"},
                )
            ]
            + _supercell,
        )

        self.dielectric_property = ipw.Dropdown(
            options=[
                ["dielectric tensor", "dielectric"],
                #'ir',
                #'raman',
                #'born-charges',
                #'nac',
                #'bec',
                #'susceptibility-derivative',
                #'non-linear-susceptibility',
                ["none", "none"],
            ],
            value="none",
            description="Dielectric property:",
            disabled=False,
            style={"description_width": "initial"},
        )

        # to trigger Dielectric property = Raman... FOR POLAR MATERIALS.
        self.material_is_polar = ipw.ToggleButtons(
            options=[("Off", "off"), ("On", "on")],
            value="off",
            style={"description_width": "initial"},
        )

        self.children = [
            self.settings_title,
            self.settings_help,
            ipw.HBox(
                children=[
                    self.phonon_property,
                    self.supercell_selector,
                ],
                layout=ipw.Layout(justify_content="flex-start"),
            ),
            self.dielectric_property,
            self.polar_help,
            ipw.HBox(
                children=[
                    ipw.Label(
                        "Material is polar:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.material_is_polar,
                ]
            ),
        ]
        super().__init__(**kwargs)

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        if isinstance(self.phonon_property, str):
            return {
                "phonon_property": self.phonon_property,
                "dielectric_property": self.dielectric_property,
                "material_is_polar": self.material_is_polar,
                "supercell_selector": self.supercell,
            }
        return {
            "phonon_property": self.phonon_property.value,
            "dielectric_property": self.dielectric_property.value,
            "material_is_polar": self.material_is_polar.value,
            "supercell_selector": self.supercell,
        }

    def load_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.phonon_property.value = input_dict.get("phonon_property", "NONE")
        self.dielectric_property.value = input_dict.get("dielectric_property", "none")
        self.material_is_polar.value = input_dict.get("material_is_polar", "off")
        self.supercell = input_dict.get("supercell_selector", [1, 1, 1])
