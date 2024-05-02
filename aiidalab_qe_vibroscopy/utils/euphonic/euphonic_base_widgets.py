import pathlib
import tempfile
import io

import base64
from IPython.display import HTML, clear_output, display

import euphonic
from phonopy.file_IO import write_force_constants_to_hdf5, write_disp_yaml

import ipywidgets as ipw
import plotly.graph_objects as go
import plotly.io as pio

# from ..euphonic.bands_pdos import *
from .intensity_maps import *
import json
from monty.json import jsanitize

# sys and os used to prevent euphonic to print in the stdout.
import sys
import os

########################
################################ START DESCRIPTION
########################

"""
In this module we have the functions and widgets to be used in the app.
Essentially we create the force constants (fc) instance via the phonopy.yaml.

def export_phononworkchain_data(node, fermi_energy=None):
Functions from intensity_maps.py and bands_pdos.py are used in order to computed the quantities, in the
export_phononworkchain_data function, used then in the result.py panel.
"""

########################
################################ END DESCRIPTION
########################

COLORSCALE = "Viridis"
COLORBAR_DICT = dict(orientation="v", showticklabels=False, x=1, thickness=10, len=0.4)

# # Intensity map widget
class StructureFactorBasePlotWidget(ipw.VBox):
    """
    Widget to plot the Structure Factor for single crystals or powder samples.
    It takes as input the spectra as generated via the
    `produce_bands_weigthed_data` or `produce_powder_data` functions, called in the
    __init__ of the master widgets, respectively: `SingleCrystalFullWidget` and `PowderFullWidget`.

    NB: Intensity is relative to the maximum intensity at T=0K.

    We use this as base for three sub widgets: intensity, powder and q-sections.
    """

    THz_to_meV = 4.13566553853599  # conversion factor.

    def __init__(self, final_xspectra, **kwargs):
        """

        Args:
            final_xspectra (_type_):
        """

        self.message_fig = ipw.HTML("")
        self.message_fig.layout.display = "none"

        if self.fig.layout.images:
            for image in self.fig.layout.images:
                image["scl"] = 2  # Set the scale for each image

        self.fig["layout"]["xaxis"].update(
            range=[min(final_xspectra), max(final_xspectra)]
        )

        self.fig["layout"]["yaxis"].update(title="meV")

        self.fig.update_layout(
            height=500,
            width=700,
            margin=dict(l=15, r=15, t=15, b=15),
        )
        # Update x-axis and y-axis to enable autoscaling
        self.fig.update_xaxes(autorange=True)
        self.fig.update_yaxes(autorange=True)

        # Update the layout to enable autoscaling
        self.fig.update_layout(autosize=True)

        self.slider_intensity = ipw.FloatRangeSlider(
            value=[1, 100],  # Default selected interval
            min=1,
            max=100,
            step=1,
            orientation="horizontal",
            readout=True,
            readout_format=".0f",
            layout=ipw.Layout(
                width="400px",
            ),
        )
        self.slider_intensity.observe(self._update_intensity_filter, "value")

        self.specification_intensity = ipw.HTML(
            "(Intensity is relative to the maximum intensity at T=0K)"
        )

        self.E_units_button = ipw.ToggleButtons(
            options=[
                ("meV", "meV"),
                ("THz", "THz"),
                # ("cm<sup>-1</sup>", "cm-1"),
            ],
            value="meV",
            description="Energy units:",
            disabled=False,
            layout=ipw.Layout(
                width="auto",
            ),
        )
        self.E_units_button.observe(self._update_energy_units, "value")

        # Create and show figure
        super().__init__(
            children=[
                self.message_fig,
                self.fig,
                ipw.HBox([ipw.HTML("Intensity window (%):"), self.slider_intensity]),
                self.specification_intensity,
                self.E_units_button,
            ],
            layout=ipw.Layout(
                width="100%",
            ),
        )

    def _update_spectra(
        self,
        final_zspectra,
    ):
        # this will be called in the _update_spectra method of SingleCrystalPlotWidget and PowderPlotWidget

        # Update the layout to enable autoscaling
        self.fig.update_layout(autosize=True)

        # We should do a check, if we have few points (<200?) provide like a warning..
        # Also decise less than what, 30%, 50%...?

        """
        visible_points = len(
            np.where(self.fig.data[0].z > 0.5)[0]
        )
        if visible_points < 1000:
            message = f"Only {visible_points}/{len(final_zspectra.T)} points have intensity higher than 50%"
            self.message_fig.value = message
            self.message_fig.layout.display = "block"
        else:
            self.message_fig.layout.display = "none"
        """

        # I have also to update the energy window. or better, to set the intensity to respect the current intensity window selected:
        self.fig.data[0].zmax = (
            self.slider_intensity.value[1] * np.max(self.fig.data[0].z) / 100
        )  # above this, it is all yellow, i.e. max intensity.
        self.fig.data[0].zmin = (
            self.slider_intensity.value[0] * np.max(self.fig.data[0].z) / 100
        )  # above this, it is all yellow, i.e. max intensity.

    def _update_intensity_filter(self, change):
        # the value of the intensity slider is in fractions of the max.
        if change["new"] != change["old"]:
            self.fig.data[0].zmax = (
                change["new"][1] * np.max(self.fig.data[0].z) / 100
            )  # above this, it is all yellow, i.e. max intensity.
            self.fig.data[0].zmin = (
                change["new"][0] * np.max(self.fig.data[0].z) / 100
            )  # below this, it is all blue, i.e. zero intensity

    def _update_energy_units(self, change):
        # the value of the intensity slider is in fractions of the max.
        if change["new"] != change["old"]:
            self.fig.data[0].y = (
                self.fig.data[0].y * self.THz_to_meV
                if change["new"] == "meV"
                else self.fig.data[0].y / self.THz_to_meV
            )

        self.fig["layout"]["yaxis"].update(title=change["new"])


#### SETTINGS WIDGET:


class StructureFactorSettingsBaseWidget(ipw.VBox):

    """
    Collects all the button and widget used to define settings for Neutron dynamic structure factor,
    both single crystal or powder.
    """

    def __init__(self, **kwargs):

        super().__init__()

        self.float_q_spacing = ipw.FloatText(
            value=0.01,
            step=0.001,
            description="q step (1/A)",
            tooltip="q spacing in 1/A",
        )
        self.float_q_spacing.observe(self._on_setting_changed, names="value")

        self.float_energy_broadening = ipw.FloatText(
            value=0.5,
            step=0.01,
            description="&Delta;E (meV)",
            tooltip="Energy broadening in meV",
        )
        self.float_energy_broadening.observe(self._on_setting_changed, names="value")

        self.int_energy_bins = ipw.IntText(
            value=200,
            description="#E bins",
            tooltip="Number of energy bins",
        )
        self.int_energy_bins.observe(self._on_setting_changed, names="value")

        self.float_T = ipw.FloatText(
            value=0,
            step=0.01,
            description="T (K)",
            disabled=False,
        )
        self.float_T.observe(self._on_setting_changed, names="value")

        self.weight_button = ipw.ToggleButtons(
            options=[
                ("Coherent", "coherent"),
                ("DOS", "dos"),
            ],
            value="coherent",
            description="weight:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.weight_button.observe(self._on_weight_button_change, names="value")

        self.plot_button = ipw.Button(
            description="Replot",
            icon="pencil",
            button_style="primary",
            disabled=True,
            layout=ipw.Layout(width="auto"),
        )
        self.plot_button.observe(self._on_plot_button_changed, names="disabled")

        self.reset_button = ipw.Button(
            description="Reset",
            icon="recycle",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(width="auto"),
        )

        self.download_button = ipw.Button(
            description="Download Data and Plot",
            icon="download",
            button_style="primary",
            disabled=False,  # Large files...
            layout=ipw.Layout(width="auto"),
        )

        self.reset_button.on_click(self._reset_settings)

    def _reset_settings(self, _):
        self.float_q_spacing.value = 0.01
        self.float_energy_broadening.value = 0.5
        self.int_energy_bins.value = 200
        self.float_T.value = 0
        self.weight_button.value = "coherent"

    def _on_plot_button_changed(self, change):
        if change["new"] != change["old"]:
            self.download_button.disabled = not change["new"]

    def _on_weight_button_change(self, change):
        if change["new"] != change["old"]:
            self.float_T.value = 0
            self.float_T.disabled = True if change["new"] == "dos" else False
            self.plot_button.disabled = False

    def _on_setting_changed(self, change):
        self.plot_button.disabled = False


class SingleCrystalSettingsWidget(StructureFactorSettingsBaseWidget):
    def __init__(self, **kwargs):

        self.custom_kpath_description = ipw.HTML(
            """
            <div style="padding-top: 0px; padding-bottom: 0px; line-height: 140%;">
                <b>Custom q-points path for the structure factor</b>: <br>
                we can provide it via a specific format: <br>
                (1) each linear path should be divided by '|'; <br>
                (2) each path is composed of 'qxi qyi qzi - qxf qyf qzf' where qxi and qxf are, respectively,
                the start and end x-coordinate of the q direction, in reciprocal lattice units (rlu).<br>
                An example path is: '0 0 0 - 1 1 1 | 1 1 1 - 0.5 0.5 0.5'. <br>
                For now, we do not support fractions (i.e. we accept 0.5 but not 1/2).
            </div>
            """
        )

        self.custom_kpath_text = ipw.Text(
            value="",
            description="Custom path (rlu):",
            style={"description_width": "initial"},
        )
        custom_style = '<style>.custom-font { font-family: "Monospaced"; font-size: 16px; }</style>'
        display(ipw.HTML(custom_style))
        self.custom_kpath_text.add_class("custom-font")

        self.custom_kpath_text.observe(self._on_setting_changed, names="value")

        # Please note: if you change the order of the widgets below, it will
        # affect the usage of the children[0] below in the full widget.

        super().__init__()

        self.children = [
            ipw.HBox(
                [
                    ipw.VBox(
                        [
                            ipw.HBox(
                                [
                                    self.reset_button,
                                    self.plot_button,
                                    self.download_button,
                                ]
                            ),
                            self.specification_intensity,
                            self.float_q_spacing,
                            self.float_energy_broadening,
                            self.int_energy_bins,
                            self.float_T,
                            self.weight_button,
                        ],
                        layout=ipw.Layout(
                            width="50%",
                        ),
                    ),
                    ipw.VBox(
                        [
                            self.custom_kpath_description,
                            self.custom_kpath_text,
                        ],
                        layout=ipw.Layout(
                            width="80%",
                        ),
                    ),
                ],  # end of HBox children
            ),
        ]

    def _reset_settings(self, _):
        self.custom_kpath_text.value = ""
        super()._reset_settings(_)
