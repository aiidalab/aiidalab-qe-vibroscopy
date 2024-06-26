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


class SingleCrystalPlotWidget(StructureFactorBasePlotWidget):
    def __init__(self, spectra, intensity_ref_0K=1, **kwargs):

        (
            final_xspectra,
            final_zspectra,
            ticks_positions,
            ticks_labels,
        ) = generated_curated_data(spectra)
        # Data to contour is the sum of two Gaussian functions.
        x, y = np.meshgrid(spectra[0].x_data.magnitude, spectra[0].y_data.magnitude)

        self.intensity_ref_0K = intensity_ref_0K

        self.fig = go.FigureWidget()

        heatmap_trace = go.Heatmap(
            z=final_zspectra.T,
            y=y[:, 0] * self.THz_to_meV,
            x=None,
            colorbar=COLORBAR_DICT,
            colorscale=COLORSCALE,  # imported from euphonic_base_widgets
        )

        # Add colorbar
        colorbar = heatmap_trace.colorbar
        colorbar.x = 1.05  # Move colorbar to the right
        colorbar.y = 0.5  # Center colorbar vertically

        # Add heatmap trace to figure
        self.fig.add_trace(heatmap_trace)

        self.fig.update_layout(
            xaxis=dict(
                tickmode="array", tickvals=ticks_positions, ticktext=ticks_labels
            )
        )
        self.fig["layout"]["yaxis"].update(range=[min(y[:, 0]), max(y[:, 0])])

        # Create and show figure
        super().__init__(
            final_xspectra,
            **kwargs,
        )

    def _update_spectra(self, spectra):

        (
            final_xspectra,
            final_zspectra,
            ticks_positions,
            ticks_labels,
        ) = generated_curated_data(spectra)
        # Data to contour is the sum of two Gaussian functions.
        x, y = np.meshgrid(spectra[0].x_data.magnitude, spectra[0].y_data.magnitude)

        # If I do this
        #   self.data = ()
        # I have a delay in the plotting, I have blank plot while it
        # is adding the new trace (see below); So, I will instead do the
        # re-assignement of the self.data = [self.data[1]] afterwards.

        x = None  # if mode == "intensity" else x[0]
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
        self.fig.update_layout(
            xaxis=dict(
                tickmode="array", tickvals=ticks_positions, ticktext=ticks_labels
            )
        )

        self.fig.data = [self.fig.data[1]]

        super()._update_spectra(final_zspectra)


class SingleCrystalSettingsWidget(StructureFactorSettingsBaseWidget):
    def __init__(self, **kwargs):

        self.custom_kpath_description = ipw.HTML(
            """
            <div style="padding-top: 0px; padding-bottom: 0px; line-height: 140%;">
                <b>Custom q-points path for the structure factor</b>: <br>
                we can provide it via a specific format: <br>
                (1) each linear path should be divided by '|'; <br>
                (2) each path is composed of 'qxi qyi qzi - qxf qyf qzf' where qxi and qxf are, respectively,
                the start and end x-coordinate of the q direction, in crystal coordinates.<br>
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
                            self.float_q_spacing,
                            self.float_energy_broadening,
                            self.int_energy_bins,
                            self.float_T,
                            self.weight_button,
                        ],
                        layout=ipw.Layout(
                            width="60%",
                        ),
                    ),
                    ipw.VBox(
                        [
                            self.custom_kpath_description,
                            self.custom_kpath_text,
                        ],
                        layout=ipw.Layout(
                            width="70%",
                        ),
                    ),
                ],  # end of HBox children
            ),
        ]

    def _reset_settings(self, _):
        self.custom_kpath_text.value = ""
        super()._reset_settings(_)


class SingleCrystalFullWidget(ipw.VBox):
    """
    The Widget to display specifically the Intensity map of Dynamical structure
    factor for single crystal.

    The scattering lengths used in the `produce_bands_weigthed_data` function
    are tabulated (Euphonic/euphonic/data/sears-1992.json)
    and are from Sears (1992) Neutron News 3(3) pp26--37.
    """

    def __init__(self, fc, **kwargs):

        self.fc = fc

        self.spectra, self.parameters = produce_bands_weigthed_data(
            fc=self.fc, plot=False  # CHANGED
        )

        self.title_intensity = ipw.HTML(
            "<h3>Neutron dynamic structure factor - Single Crystal</h3>"
        )

        self.settings_intensity = SingleCrystalSettingsWidget()
        self.settings_intensity.plot_button.on_click(self._on_plot_button_clicked)
        self.settings_intensity.download_button.on_click(self.download_data)

        # This is used in order to have an overall intensity scale.
        self.intensity_ref_0K = np.max(self.spectra[0].z_data.magnitude)  # CHANGED

        self.map_widget = SingleCrystalPlotWidget(
            self.spectra, intensity_ref_0K=self.intensity_ref_0K
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
            }
        )
        parameters_ = AttrDict(self.parameters)  # CHANGED

        # custom linear path
        if len(self.settings_intensity.custom_kpath_text.value) > 1:
            coordinates, labels = self.curate_path_and_labels()
            linear_path = {
                "coordinates": coordinates,
                "labels": labels,  # ["$\Gamma$","X","X","(1,1,1)"],
                "delta_q": parameters_["q_spacing"],
            }
        else:
            linear_path = None

        self.spectra, self.parameters = produce_bands_weigthed_data(
            params=parameters_,
            fc=self.fc,
            plot=False,
            linear_path=linear_path,  # CHANGED
        )

        self.settings_intensity.plot_button.disabled = True
        self.map_widget._update_spectra(self.spectra)  # CHANGED

    def download_data(self, _=None):
        """
        Download both the ForceConstants and the spectra json files.
        """
        force_constants_dict = self.fc.to_dict()

        filename = "single_crystal.json"
        my_dict = {}
        for branch in range(len(self.spectra)):
            my_dict[str(branch)] = self.spectra[branch].to_dict()
        my_dict.update(
            {
                "weighting": self.settings_intensity.weight_button.value,
                "q_spacing": self.settings_intensity.float_q_spacing.value,
                "energy_broadening": self.settings_intensity.float_energy_broadening.value,
                "ebins": self.settings_intensity.int_energy_bins.value,
                "temperature": self.settings_intensity.float_T.value,
            }
        )
        for k in ["weighting", "q_spacing", "temperature"]:
            filename += "_" + k + "_" + str(my_dict[k])

        # FC download:
        json_str = json.dumps(jsanitize(force_constants_dict))
        b64_str = base64.b64encode(json_str.encode()).decode()
        self._download(payload=b64_str, filename="force_constants.json")

        # Powder data download:
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

    def curate_path_and_labels(
        self,
    ):
        # I do not like this implementation (MB)
        coordinates = []
        labels = []
        path = self.settings_intensity.custom_kpath_text.value
        linear_paths = path.split("|")
        for i in linear_paths:
            scoords = []
            s = i.split(
                " - "
            )  # not i.split("-"), otherwise also the minus of the negative numbers are used for the splitting.
            for k in s:
                labels.append(k.strip())
                # AAA missing support for fractions.
                l = tuple(map(float, [kk for kk in k.strip().split(" ")]))
                scoords.append(l)
            coordinates.append(scoords)
        return coordinates, labels

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
