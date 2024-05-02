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
from ..euphonic.intensity_maps import *

import json
from monty.json import jsanitize

# sys and os used to prevent euphonic to print in the stdout.
import sys
import os

from aiidalab_qe_vibroscopy.utils.euphonic.euphonic_base_widgets import *


class PowderPlotWidget(StructureFactorBasePlotWidget):
    def __init__(self, spectra, intensity_ref_0K=1, **kwargs):

        final_zspectra = spectra.z_data.magnitude
        final_xspectra = spectra.x_data.magnitude
        # Data to contour is the sum of two Gaussian functions.
        x, y = np.meshgrid(spectra.x_data.magnitude, spectra.y_data.magnitude)

        self.intensity_ref_0K = intensity_ref_0K

        self.fig = go.FigureWidget()

        self.fig.add_trace(
            go.Heatmap(
                z=final_zspectra.T,
                y=y[:, 0] * self.THz_to_meV,
                x=x[0],
                colorbar=COLORBAR_DICT,
                colorscale=COLORSCALE,  # imported from euphonic_base_widgets
            )
        )

        self.fig["layout"]["xaxis"].update(title="|q| (1/A)")
        self.fig["layout"]["yaxis"].update(range=[min(y[:, 0]), max(y[:, 0])])

        # Create and show figure
        super().__init__(
            final_xspectra,
            **kwargs,
        )

    def _update_spectra(self, spectra):

        final_zspectra = spectra.z_data.magnitude
        final_xspectra = spectra.x_data.magnitude
        # Data to contour is the sum of two Gaussian functions.
        x, y = np.meshgrid(spectra.x_data.magnitude, spectra.y_data.magnitude)

        # If I do this
        #   self.data = ()
        # I have a delay in the plotting, I have blank plot while it
        # is adding the new trace (see below); So, I will instead do the
        # re-assignement of the self.data = [self.data[1]] afterwards.

        x = x[0]
        self.fig.add_trace(
            go.Heatmap(
                z=final_zspectra.T,
                y=y[:, 0] * self.THz_to_meV
                if self.E_units_button.value == "meV"
                else y[:, 0],
                x=x,
                colorbar=COLORBAR_DICT,
                colorscale=COLORSCALE,  # imported from euphonic_base_widgets
            )
        )

        # change the path wants also a change in the labels.
        # this is delays things
        self.fig["layout"]["xaxis"].update(title="|q| (1/A)")

        self.fig.data = [self.fig.data[1]]

        super()._update_spectra(final_zspectra)


class PowderSettingsWidget(StructureFactorSettingsBaseWidget):
    def __init__(self, **kwargs):

        self.float_qmin = ipw.FloatText(
            value=0,
            description="|q|<sub>min</sub> (1/A)",
        )
        self.float_qmin.observe(self._on_setting_changed, names="value")

        self.float_qmax = ipw.FloatText(
            step=0.01,
            value=1,
            description="|q|<sub>max</sub> (1/A)",
        )
        self.float_qmax.observe(self._on_setting_changed, names="value")

        self.int_npts = ipw.IntText(
            value=100,
            description="npts",
            tooltip="Number of points to be used in the average sphere.",
        )
        self.int_npts.observe(self._on_setting_changed, names="value")

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
                            self.float_q_spacing,
                            self.float_qmin,
                            self.float_qmax,
                            # self.int_npts,
                        ],
                        layout=ipw.Layout(
                            width="50%",
                        ),
                    ),
                    ipw.VBox(
                        [
                            self.float_energy_broadening,
                            self.int_energy_bins,
                            self.float_T,
                            self.weight_button,
                        ],
                        layout=ipw.Layout(
                            width="50%",
                        ),
                    ),
                ],  # end of HBox children
            ),
        ]

    def _reset_settings(self, _):
        self.float_qmin.value = 0
        self.float_qmax.value = 1
        self.int_npts.value = 100
        super()._reset_settings(_)


class PowderFullWidget(ipw.VBox):
    """
    The Widget to display specifically the Intensity map of Dynamical structure
    factor for powder samples.

    The scattering lengths used in the `produce_bands_weigthed_data` function
    are tabulated (Euphonic/euphonic/data/sears-1992.json)
    and are from Sears (1992) Neutron News 3(3) pp26--37.
    """

    def __init__(self, fc, intensity_ref_0K=1, **kwargs):

        self.fc = fc

        self.spectra, self.parameters = produce_powder_data(
            params=parameters_powder, fc=self.fc
        )

        self.title_intensity = ipw.HTML(
            "<h3>Neutron dynamic structure factor - Powder sample</h3>"
        )

        self.settings_intensity = PowderSettingsWidget(mode="powder")
        self.settings_intensity.plot_button.on_click(self._on_plot_button_clicked)
        self.settings_intensity.download_button.on_click(self.download_data)

        # This is used in order to have an overall intensity scale. Inherithed from the SingleCrystal
        self.intensity_ref_0K = intensity_ref_0K  # CHANGED

        self.map_widget = PowderPlotWidget(
            self.spectra, mode="powder", intensity_ref_0K=self.intensity_ref_0K
        )  # CHANGED

        super().__init__(
            children=[
                self.title_intensity,
                self.map_widget,
                self.settings_intensity,
            ],
        )

    def _on_plot_button_clicked(self, change=None):
        self.parameters.update(
            {
                "weighting": self.settings_intensity.weight_button.value,
                "q_spacing": self.settings_intensity.float_q_spacing.value,
                "energy_broadening": self.settings_intensity.float_energy_broadening.value,
                "ebins": self.settings_intensity.int_energy_bins.value,
                "temperature": self.settings_intensity.float_T.value,
                "q_min": self.settings_intensity.float_qmin.value,
                "q_max": self.settings_intensity.float_qmax.value,
                "npts": self.settings_intensity.int_npts.value,
            }
        )
        parameters_powder = AttrDict(self.parameters)

        self.spectra, self.parameters = produce_powder_data(
            params=parameters_powder, fc=self.fc, plot=False
        )

        self.settings_intensity.plot_button.disabled = True
        self.map_widget._update_spectra(self.spectra)  # CHANGED

    def download_data(self, _=None):
        """
        Download both the ForceConstants and the spectra json files.
        """
        force_constants_dict = self.fc.to_dict()

        filename = "powder"
        my_dict = self.spectra.to_dict()
        my_dict.update(
            {
                "weighting": self.settings_intensity.weight_button.value,
                "q_spacing": self.settings_intensity.float_q_spacing.value,
                "energy_broadening": self.settings_intensity.float_energy_broadening.value,
                "ebins": self.settings_intensity.int_energy_bins.value,
                "temperature": self.settings_intensity.float_T.value,
                "q_min": self.settings_intensity.float_qmin.value,
                "q_max": self.settings_intensity.float_qmax.value,
                # "npts": self.settings_intensity.int_npts.value,
            }
        )
        for k in ["weighting", "q_spacing", "temperature"]:
            filename += "_" + k + "_" + str(my_dict[k])

        # FC download:
        json_str = json.dumps(jsanitize(force_constants_dict))
        b64_str = base64.b64encode(json_str.encode()).decode()
        self._download(payload=b64_str, filename="force_constants.json")

        # Powder download:
        json_str = json.dumps(jsanitize(my_dict))
        b64_str = base64.b64encode(json_str.encode()).decode()
        self._download(payload=b64_str, filename=filename + ".json")

        # Plot download:
        ## Convert the FigureWidget to an image in base64 format
        image_bytes = pio.to_image(
            self.map_widget.children[1], format="png", width=800, height=600
        )
        b64_str = base64.b64encode(image_bytes).decode()
        self._download(payload=b64_str, filename=filename + ".png")

    @staticmethod
    def _download(payload, filename):
        from IPython.display import Javascript

        javas = Javascript(
            """
            var link = document.createElement('a');
            link.href = 'data:text/json;charset=utf-8;base64,{payload}'
            link.download = "{filename}"
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            """.format(
                payload=payload, filename=filename
            )
        )
        display(javas)
