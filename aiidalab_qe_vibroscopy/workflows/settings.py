# -*- coding: utf-8 -*-
"""Panel for PhononWorkchain plugin.

Authors:

    * Miki Bonacci <miki.bonacci@psi.ch>
    Inspired by Xing Wang <xing.wang@psi.ch>
"""
import ipywidgets as ipw
from aiida.orm import Float, Int, Str
import traitlets as tl

from aiidalab_qe.common.panel import Panel


class Setting(Panel):
    title = "Vibrational Settings"

    def __init__(self, **kwargs):
        self.settings_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h4>Phonons and dielectric settings</h4></div>"""
        )
        self.settings_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Please select the phononic and dielectric properties to be computed in the simulation. Usually,
            supercell size of 2x2x2 is at convergence with respect to finite atomic displacement in the unit cell.
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

        # I want to be able to select more than only one... this has to change at the PhononWorkChain level.
        self.phonon_property_ = ipw.Dropdown(
            options=[
                ["band structure", "BANDS"],
                ["density of states (DOS)", "DOS"],
                ["thermal properties", "THERMODYNAMIC"],
                ["force constants", "NONE"],
                ["none", "none"],
            ],
            value="none",
            description="Phonon property:",
            disabled=False,
            style={"description_width": "initial"},
        )

        self.dielectric_property_ = ipw.Dropdown(
            options=[
                ["dielectric tensor", "dielectric"],
                ["infrared", "ir"],
                ["raman", "raman"],
                ["born-charges", "born-charges"],
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

        self.spectrum_ = ipw.ToggleButtons(
            options=[("Off", "off"), ("Infrared", "ir"), ("Raman", "raman")],
            value="off",
            style={"description_width": "initial"},
        )

        self.spectrum_.observe(self._onclick_spectrum_choice, "value")

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
        self.material_is_polar_.observe(self._onclick_material_is_polar, "value")

        self.children = [
            self.settings_title,
            self.settings_help,
            ipw.HBox(
                children=[
                    ipw.Label(
                        "Spectroscopy:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.spectrum_,
                ]
            ),
            ipw.HBox(
                children=[
                    self.phonon_property_,
                    self.supercell_widget,
                ],
                layout=ipw.Layout(justify_content="flex-start"),
            ),
            self.dielectric_property_,
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

    def _onclick_spectrum_choice(self, change):
        # if spectra are choosen, we change the dielectric property
        # and we do simulation on 1x1x1.
        if change["new"] in ["raman", "ir"]:
            self.dielectric_property_.value = change["new"]
            self.phonon_property_.value = "none"
            self.dielectric_property_.disabled = True
            self.phonon_property_.disabled = True
            for elem in [self._sc_x, self._sc_y, self._sc_z]:
                elem.value = 1
                elem.disabled = True
        elif change["new"] == "off":
            self.dielectric_property_.disabled = False
            self.phonon_property_.disabled = False

            # logic for polar materials, if we turn off spectra:
            if self.material_is_polar_.value == "on":
                self.dielectric_property_.value = "raman"
            else:
                self.dielectric_property_.value = "none"

            for elem in [self._sc_x, self._sc_y, self._sc_z]:
                elem.disabled = False

    def _onclick_material_is_polar(self, change):
        # if spectra are choosen, we change the dielectric property
        # and we do simulation on 1x1x1.
        if change["new"] == "on" and self.spectrum_.value == "off":
            self.dielectric_property_.value = "raman"
            self.dielectric_property_.disabled = True
        elif change["new"] == "off" and self.spectrum_.value == "off":
            self.dielectric_property_.disabled = False

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        return {
            "phonon_property": self.phonon_property_.value,
            "dielectric_property": self.dielectric_property_.value,
            "material_is_polar": self.material_is_polar_.value,
            "supercell_selector": self.supercell,
            "spectrum": self.spectrum_.value,
        }

    def load_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.phonon_property_.value = input_dict.get("phonon_property", "none")
        self.dielectric_property_.value = input_dict.get("dielectric_property", "none")
        self.material_is_polar_.value = input_dict.get("material_is_polar", "off")
        self.spectrum_.value = input_dict.get("spectrum", "off")
        self.supercell = input_dict.get("supercell_selector", [2, 2, 2])

    def reset(self):
        """Reset the panel"""
        self.phonon_property_.value = "none"
        self.dielectric_property_.value = "none"
        self.material_is_polar_.value = "off"
        self.spectrum_.value = "off"
        self.supercell = [2, 2, 2]
