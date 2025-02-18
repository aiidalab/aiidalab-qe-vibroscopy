# -*- coding: utf-8 -*-
"""Panel for PhononWorkchain plugin.

Authors:

    * Miki Bonacci <miki.bonacci@psi.ch>
    Inspired by Xing Wang <xing.wang@psi.ch>
"""

import ipywidgets as ipw

from aiidalab_qe.common.panel import ConfigurationSettingsPanel
from aiidalab_qe_vibroscopy.app.model import VibroConfigurationSettingsModel

from aiida.plugins import DataFactory

HubbardStructureData = DataFactory("quantumespresso.hubbard_structure")


class VibroConfigurationSettingPanel(
    ConfigurationSettingsPanel[VibroConfigurationSettingsModel],
):
    title = "Vibrational Settings"
    identifier = "vibronic"

    def __init__(self, model: VibroConfigurationSettingsModel, **kwargs):
        super().__init__(model, **kwargs)

        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )

    def render(self):
        if self.rendered:
            return

        self.settings_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h4>Vibrational Settings</h4></div>"""
        )
        self.settings_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            <ul>
                <li>Calculations are performed using the <b><a href="https://aiida-vibroscopy.readthedocs.io/en/latest/"
                target="_blank">aiida-vibroscopy</b></a> plugin (L. Bastonero and N. Marzari, <a href="https://www.nature.com/articles/s41524-024-01236-3"
                target="_blank">npj Comput. Mater. <b>10</b>, 55, 2024</a>).</li>
                <ul>
                    <li>The plugin employes the finite-displacement and finite-field approach.</li>
                    <li>Raman spectra are simulated in the first-order non-resonant regime. </li>
                </ul>
            <li>The inelastic neutron scattering
            structure factor is calculated as post processing using the <b><a href="https://euphonic.readthedocs.io/en/stable/index.html#" target="_blank">Euphonic</b></a> code (R. Fair et al.,
            <a href="https://doi.org/10.1107/S1600576722009256"
            target="_blank">J. Appl. Cryst. <b>55</b>, 1689, 2022</a>).</li>
            </ul>
            </div>""",
            layout=ipw.Layout(width="400"),
        )

        self.use_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h5>Available simulations:</h5></div>"""
        )

        self.use_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            <ul>
                <li> <em>IR/Raman spectra</em>: both single crystal and powder samples.</li>
                <li> <em>Phonons properties</em>: bands, density
                of states and thermal properties (Helmoltz free energy, entropy and specific heat at constant volume).</li>
                <li> <em>Dielectric properties</em>: Born charges,
                high-frequency dielectric tensor, non-linear optical susceptibility and raman tensors.</li>
                <li> <em>Inelastic neutron scattering (INS)</em>: dynamic structure factor and powder intensity maps.</li>
            </ul>
            </div>""",
            layout=ipw.Layout(width="400"),
        )

        self.hint_button_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Select a supercell size for Phonon properties:
            <ul>
                <li>Larger supercells increase computational costs.</li>
                <li>A 2x2x2 supercell is usually adequate.</li>
            </ul>
            You can use the <em>Size hint</em> button for an estimate, performed imposing a minimum lattice vector magnitude of 15Å along the periodic directions.
            </div>""",
        )

        self.simulation_type = ipw.Dropdown(
            layout=ipw.Layout(width="450px"),
        )
        ipw.dlink(
            (self._model, "simulation_type_options"),
            (self.simulation_type, "options"),
        )
        ipw.link(
            (self._model, "simulation_type"),
            (self.simulation_type, "value"),
        )
        self.simulation_type.observe(
            self._on_change_simulation_type,
            "value",
        )

        self.symmetry_symprec = ipw.BoundedFloatText(
            max=1,
            min=1e-7,
            step=1e-4,
            description="Symmetry tolerance (symprec):",
            style={"description_width": "initial"},
            layout={"width": "300px"},
        )
        ipw.link(
            (self._model, "symmetry_symprec"),
            (self.symmetry_symprec, "value"),
        )

        self.supercell_x = ipw.BoundedIntText(
            min=1,
            layout={"width": "40px"},
        )
        self.supercell_y = ipw.BoundedIntText(
            min=1,
            layout={"width": "40px"},
        )
        self.supercell_z = ipw.BoundedIntText(
            min=1,
            layout={"width": "40px"},
        )
        ipw.link(
            (self._model, "supercell_x"),
            (self.supercell_x, "value"),
        )
        ipw.link(
            (self._model, "disable_x"),
            (self.supercell_x, "disabled"),
        )
        ipw.link(
            (self._model, "supercell_y"),
            (self.supercell_y, "value"),
        )
        ipw.link(
            (self._model, "disable_y"),
            (self.supercell_y, "disabled"),
        )
        ipw.link(
            (self._model, "supercell_z"),
            (self.supercell_z, "value"),
        )
        ipw.link(
            (self._model, "disable_z"),
            (self.supercell_z, "disabled"),
        )

        self.supercell_selector = ipw.HBox(
            children=[
                ipw.HTML(
                    description="Supercell size:",
                    style={"description_width": "initial"},
                )
            ]
            + [
                self.supercell_x,
                self.supercell_y,
                self.supercell_z,
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
        self.supercell_hint_button.on_click(self._model.suggest_supercell)

        # reset supercell
        self.supercell_reset_button = ipw.Button(
            description="Reset hint",
            disabled=False,
            layout=ipw.Layout(width="100px"),
            button_style="warning",
        )
        # supercell reset reaction
        self.supercell_reset_button.on_click(self._model.supercell_reset)

        self.symmetry_symprec_reset_button = ipw.Button(
            description="Reset symprec",
            disabled=False,
            layout=ipw.Layout(width="125px"),
            button_style="warning",
        )
        # supercell reset reaction
        self.symmetry_symprec_reset_button.on_click(self._model.reset_symprec)

        # Estimate the number of supercells for frozen phonons.
        self.supercell_number_estimator = ipw.HTML(
            style={"description_width": "initial"},
            layout=ipw.Layout(display="none"),
        )
        ipw.link(
            (self._model, "supercell_number_estimator"),
            (self.supercell_number_estimator, "value"),
        )

        # Estimate supercell button
        self.supercell_estimate_button = ipw.Button(
            description="Estimate number of supercells ➡",
            disabled=False,
            layout=ipw.Layout(width="240px", display="none"),
            button_style="info",
            tooltip="Number of supercells for phonons calculations;\nwarning: for large systems, this may take some time.",
        )

        # supercell reset reaction
        self.supercell_estimate_button.on_click(self._model._estimate_supercells)

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

        # self.symmetry_symprec.observe(self._activate_estimate_supercells, "value")

        # reset supercell

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
                    ipw.HTML("Select calculation:"),
                    self.simulation_type,
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

        self.rendered = True
        self._on_change_simulation_type({"new": 1})

    def _on_input_structure_change(self, _):
        self.refresh(specific="structure")
        self._model.on_input_structure_change()

    def _on_change_simulation_type(self, _):
        self.supercell_widget.layout.display = (
            "block" if self._model.simulation_type in [1, 3] else "none"
        )
        self.supercell_number_estimator.layout.display = (
            "block"
            if self._model.simulation_type in [1, 3]
            and len(self._model.input_structure.sites) <= 30
            else "none"
        )
        self.supercell_estimate_button.layout.display = (
            "block"
            if self._model.simulation_type in [1, 3]
            and len(self._model.input_structure.sites) <= 30
            else "none"
        )
