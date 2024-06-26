# -*- coding: utf-8 -*-
"""Panel for PhononWorkchain plugin.

Authors:

    * Miki Bonacci <miki.bonacci@psi.ch>
    Inspired by Xing Wang <xing.wang@psi.ch>
"""
import ipywidgets as ipw
from aiida.orm import Float, Int, Str

from aiidalab_qe.common.panel import Panel
from aiida_vibroscopy.common.properties import PhononProperty


class Setting(Panel):
    title = "Vibrational Settings"

    def __init__(self, **kwargs):
        self.settings_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h4>Vibrational spectra settings</h4></div>"""
        )
        self.settings_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Please select the vibrational spectrum to be computed in this simulation.
            </div>"""
        )
        self.workchain_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
            value="moderate",
        )

        self.spectrum = ipw.ToggleButtons(
            options=[("Infrared", "ir"), ("Raman", "raman")],
            value="raman",
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

        self.children = [
            self.settings_title,
            self.settings_help,
            ipw.HBox(
                children=[
                    ipw.Label(
                        "Spectroscopy:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.spectrum,
                ]
            ),
            self.supercell_selector,
        ]
        super().__init__(**kwargs)

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        return {"spectrum": self.spectrum.value, "supercell_selector": self.supercell}

    def load_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.spectrum.value = input_dict.get("spectrum", "raman")
        self.supercell = input_dict.get("supercell_selector", [1, 1, 1])
