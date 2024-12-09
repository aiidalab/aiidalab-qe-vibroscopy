import numpy as np
import ipywidgets as ipw
import plotly.graph_objects as go

# from ..euphonic.bands_pdos import *
from aiidalab_qe_vibroscopy.utils.euphonic.data.structure_factors import *  # noqa: F403

# sys and os used to prevent euphonic to print in the stdout.

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
    THz_to_cm1 = 33.3564095198155  # conversion factor.

    def __init__(self, model):
        super().__init__()
        self._model = model
        self.rendered = False

    def render(self):
        """Render the widget.
        This is the generic render method which can be overwritten by the subwidgets.
        However, it is important to call this method at the start of the subwidgets.render() in order to have the go.FigureWidget.
        """

        if self.rendered:
            return

        if not hasattr(self._model, "fc"):
            self._model.fetch_data()
        self._model._update_spectra()

        self.fig = go.FigureWidget()

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
        # MAYBE WE LINK ALSO THIS TO THE MODEL? so we can download the data with the preferred units.

        # Create and show figure
        self.children = [
            self.fig,
            ipw.HBox([ipw.HTML("Intensity window (%):"), self.slider_intensity]),
            self.specification_intensity,
            self.E_units_button,
        ]

    def _update_plot(self):
        """This is the generic update_plot method which can be overwritten by the subwidgets.
        However, it is important to call this method at the end of the subwidgets._update_plot() in order to update the intensity window.
        """

        self.fig.update_layout(autosize=True)

        # I have also to update the energy window. or better, to set the intensity to respect the current intensity window selected:
        self.fig.data[0].zmax = (
            self.slider_intensity.value[1] * np.max(self.fig.data[0].z) / 100
        )  # above this, it is all yellow, i.e. this is the max detachable intensity.
        self.fig.data[0].zmin = (
            self.slider_intensity.value[0] * np.max(self.fig.data[0].z) / 100
        )  # below this, it is all dark blue, i.e. this is the min detachable intensity.

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

    def __init__(self, model, **kwargs):
        super().__init__()
        self._model = model
        self.rendered = False

    def render(self):
        """Render the widget."""

        if self.rendered:
            return

        self.q_spacing = ipw.FloatText(
            value=self._model.q_spacing,
            step=0.001,
            description="q step (1/A)",
            tooltip="q spacing in 1/A",
        )
        ipw.link(
            (self._model, "q_spacing"),
            (self.q_spacing, "value"),
        )
        self.q_spacing.observe(self._on_setting_changed, names="value")

        self.energy_broadening = ipw.FloatText(
            value=self._model.energy_broadening,
            step=0.01,
            description="&Delta;E (meV)",
            tooltip="Energy broadening in meV",
        )
        ipw.link(
            (self._model, "energy_broadening"),
            (self.energy_broadening, "value"),
        )
        self.energy_broadening.observe(self._on_setting_changed, names="value")

        self.energy_bins = ipw.IntText(
            value=self._model.energy_bins,
            description="#E bins",
            tooltip="Number of energy bins",
        )
        ipw.link(
            (self._model, "energy_bins"),
            (self.energy_bins, "value"),
        )
        self.energy_bins.observe(self._on_setting_changed, names="value")

        self.temperature = ipw.FloatText(
            value=self._model.temperature,
            step=0.01,
            description="T (K)",
            disabled=False,
        )
        ipw.link(
            (self._model, "temperature"),
            (self.temperature, "value"),
        )
        self.temperature.observe(self._on_setting_changed, names="value")

        self.weight_button = ipw.ToggleButtons(
            options=[
                ("Coherent", "coherent"),
                ("DOS", "dos"),
            ],
            value=self._model.weighting,
            description="weight:",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "weighting"),
            (self.weight_button, "value"),
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
        self.reset_button.on_click(self._reset_settings)

        self.download_button = ipw.Button(
            description="Download Data and Plot",
            icon="download",
            button_style="primary",
            disabled=False,  # Large files...
            layout=ipw.Layout(width="auto"),
        )

    def _on_plot_button_changed(self, change):
        if change["new"] != change["old"]:
            self.download_button.disabled = not change["new"]

    def _on_weight_button_change(self, change):
        if change["new"] != change["old"]:
            self._model.temperature = 0
            self.temperature.disabled = True if change["new"] == "dos" else False
            self.plot_button.disabled = False

    def _on_setting_changed(
        self, change
    ):  # think if we want to do something more evident...
        self.plot_button.disabled = False

    def _reset_settings(self, _):
        self._model.reset()
