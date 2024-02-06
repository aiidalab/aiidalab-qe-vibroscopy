import pathlib
import tempfile

import base64
from IPython.display import HTML, clear_output, display

import euphonic
from phonopy.file_IO import write_force_constants_to_hdf5, write_disp_yaml

import ipywidgets as ipw
import plotly.graph_objects as go

from ..euphonic.bands_pdos import *
from ..euphonic.intensity_maps import *

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


def generate_force_constant_instance(phonopy_calc):
    blockPrint()
    ####### This is almost copied from PhonopyCalculation
    from phonopy.interface.phonopy_yaml import PhonopyYaml

    kwargs = {}

    if "settings" in phonopy_calc.inputs:
        the_settings = phonopy_calc.inputs.settings.get_dict()
        for key in ["symmetrize_nac", "factor_nac", "subtract_residual_forces"]:
            if key in the_settings:
                kwargs.update({key: the_settings[key]})

    if "phonopy_data" in phonopy_calc.inputs:
        ph = phonopy_calc.inputs.phonopy_data.get_phonopy_instance(**kwargs)
        p2s_map = phonopy_calc.inputs.phonopy_data.get_cells_mappings()["primitive"][
            "p2s_map"
        ]
        ph.produce_force_constants()
    elif "force_constants" in phonopy_calc.inputs:
        ph = phonopy_calc.inputs.force_constants.get_phonopy_instance(**kwargs)
        p2s_map = phonopy_calc.inputs.force_constants.get_cells_mappings()["primitive"][
            "p2s_map"
        ]
        ph.force_constants = phonopy_calc.inputs.force_constants.get_array(
            "force_constants"
        )

    #######

    # Create temporary directory
    #
    with tempfile.TemporaryDirectory() as dirpath:
        # phonopy.yaml generation:
        phpy_yaml = PhonopyYaml()
        phpy_yaml.set_phonon_info(ph)
        phpy_yaml_txt = str(phpy_yaml)

        with open(
            pathlib.Path(dirpath) / "phonopy.yaml", "w", encoding="utf8"
        ) as handle:
            handle.write(phpy_yaml_txt)

        # Force constants hdf5 file generation:
        # all this is needed to load the euphonic instance, in case no FC are written in phonopy.yaml

        write_force_constants_to_hdf5(
            force_constants=ph.force_constants,
            filename=pathlib.Path(dirpath) / "fc.hdf5",
            p2s_map=p2s_map,
        )

        # Read force constants (fc.hdf5) and summary+NAC (phonopy.yaml)

        fc = euphonic.ForceConstants.from_phonopy(
            path=dirpath,
            summary_name="phonopy.yaml",
            fc_name="fc.hdf5",
        )
        # print(filename)
        # print(dirpath)
    enablePrint()
    return fc


def export_euphonic_data(node, fermi_energy=None):

    if not "vibronic" in node.outputs:
        # Not a phonon calculation
        return None
    else:
        if not "phonon_bands" in node.outputs.vibronic:
            return None

    output_set = node.outputs.vibronic.phonon_bands

    phonopy_calc = output_set.creator
    fc = generate_force_constant_instance(phonopy_calc)
    # bands = compute_bands(fc)
    # pdos = compute_pdos(fc)
    return {
        "fc": fc,
    }  # "bands": bands, "pdos": pdos, "thermal": None}


def generated_curated_data(spectra):
    # here we concatenate the bands groups and create the ticks and labels.

    ticks_positions = []
    ticks_labels = []

    final_xspectra = spectra[0].x_data.magnitude
    final_zspectra = spectra[0].z_data.magnitude
    for i in spectra[1:]:
        final_xspectra = np.concatenate((final_xspectra, i.x_data.magnitude), axis=0)
        final_zspectra = np.concatenate((final_zspectra, i.z_data.magnitude), axis=0)

    for j in spectra[:]:

        # each spectra has the .x_tick_labels attribute, for the phonon bands.
        for k in j.x_tick_labels:
            ticks_positions.append(k[0])
            # ticks_labels.append("Gamma") if k[1] == '$\\Gamma$' else ticks_labels.append(k[1])
            ticks_labels.append(k[1])

            # Here below we check if we are strarting a new group, i.e. if the xticks count is starting again from 0
            if len(ticks_positions) > 3:
                if ticks_positions[-1] < ticks_positions[-2]:
                    if ticks_positions[-1] == 0:
                        ticks_positions.pop()
                        last = ticks_labels.pop()
                        ticks_labels[-1] = ticks_labels[-1] + "|" + last

                    else:
                        ticks_positions[-1] = ticks_positions[-1] + ticks_positions[-2]

    return final_xspectra, final_zspectra, ticks_positions, ticks_labels


class IntensityMapWidget(go.FigureWidget):
    """
    It is used also for powder maps.
    """

    def __init__(self, spectra, mode="intensity", **kwargs):

        # Create and show figure
        super().__init__()

        if mode == "intensity":
            (
                final_xspectra,
                final_zspectra,
                ticks_positions,
                ticks_labels,
            ) = generated_curated_data(spectra)
            # Data to contour is the sum of two Gaussian functions.
            x, y = np.meshgrid(spectra[0].x_data.magnitude, spectra[0].y_data.magnitude)
        else:
            final_zspectra = spectra.z_data.magnitude
            final_xspectra = spectra.x_data.magnitude
            # Data to contour is the sum of two Gaussian functions.
            x, y = np.meshgrid(spectra.x_data.magnitude, spectra.y_data.magnitude)

        self.add_trace(
            go.Heatmap(
                z=final_zspectra.T,
                colorbar=dict(
                    title="arb.units",
                    titleside="top",
                    tickmode="array",
                    tickvals=[2, 50, 100],
                    ticktext=["Cool", "Mild", "Hot"],
                    ticks="outside",
                ),
                y=y[:, 0],
                x=None if mode == "intensity" else x[0],
            )
        )

        if mode == "intensity":
            self.update_layout(
                xaxis=dict(
                    tickmode="array", tickvals=ticks_positions, ticktext=ticks_labels
                )
            )
        else:
            self["layout"]["xaxis"].update(title="|q| (1/A)")

        self["layout"]["xaxis"].update(range=[min(final_xspectra), max(final_xspectra)])
        self["layout"]["yaxis"].update(range=[min(y[:, 0]), max(y[:, 0])])
        self["layout"]["yaxis"].update(title="meV")

        self.update_layout(autosize=True)

    def _update_spectra(self, spectra, mode="intensity"):
        if mode == "intensity":
            (
                final_xspectra,
                final_zspectra,
                ticks_positions,
                ticks_labels,
            ) = generated_curated_data(spectra)
            # Data to contour is the sum of two Gaussian functions.
            x, y = np.meshgrid(spectra[0].x_data.magnitude, spectra[0].y_data.magnitude)
        else:
            final_zspectra = spectra.z_data.magnitude
            final_xspectra = spectra.x_data.magnitude
            # Data to contour is the sum of two Gaussian functions.
            x, y = np.meshgrid(spectra.x_data.magnitude, spectra.y_data.magnitude)

        self.data = ()
        self.add_trace(
            go.Heatmap(
                z=final_zspectra.T,
                colorbar=dict(
                    title="arb. units",
                    titleside="top",
                    tickmode="array",
                    tickvals=[2, 50, 100],
                    ticktext=["Cool", "Mild", "Hot"],
                    ticks="outside",
                ),
                y=y[:, 0],
                x=None if mode == "intensity" else x[0],
            )
        )


class IntensityFullWidget(ipw.HBox):
    def __init__(self, fc, **kwargs):

        self.fc = fc

        self.spectra, self.parameters = produce_bands_weigthed_data(
            params=parameters, fc=self.fc, plot=False
        )

        self.settings_intensity = IntensitySettingsWidget()
        self.settings_intensity.plot_button.on_click(self._on_plot_button_clicked)
        self.settings_intensity.download_button.on_click(self.download_data)

        self.fig = IntensityMapWidget(self.spectra)

        super().__init__(
            children=[
                self.settings_intensity,
                self.fig,
            ],
        )

    def _on_plot_button_clicked(self, change=None):
        self.parameters.update(
            {
                "weighting": self.children[0].weight_button.value,
                "q_spacing": self.children[0].slider_q_spacing.value,
                "energy_broadening": self.children[0].slider_energy_broadening.value,
                "ebins": self.children[0].slider_energy_bins.value,
                "temperature": self.children[0].slider_T.value,
            }
        )
        parameters = AttrDict(par_dict)

        self.spectra, self.parameters = produce_bands_weigthed_data(
            params=parameters, fc=self.fc, plot=False
        )
        self.fig._update_spectra(self.spectra)

    def download_data(self, _=None):
        filename = "intensity_spectra.json"
        my_dict = {}
        for branch in range(len(self.spectra)):
            my_dict[str(branch)] = self.spectra[branch].to_dict()
        my_dict.update(
            {
                "weighting": self.children[0].weight_button.value,
                "q_spacing": self.children[0].slider_q_spacing.value,
                "energy_broadening": self.children[0].slider_energy_broadening.value,
                "ebins": self.children[0].slider_energy_bins.value,
                "temperature": self.children[0].slider_T.value,
            }
        )
        json_str = json.dumps(jsanitize(my_dict))
        b64_str = base64.b64encode(json_str.encode()).decode()
        self._download(payload=b64_str, filename=filename)

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


class PowderFullWidget(ipw.HBox):
    def __init__(self, fc, **kwargs):

        self.fc = fc

        self.spectra, self.parameters = produce_powder_data(
            params=parameters, fc=self.fc
        )

        self.settings_intensity = IntensitySettingsWidget(mode="powder")
        self.settings_intensity.plot_button.on_click(self._on_plot_button_clicked)
        self.settings_intensity.download_button.on_click(self.download_data)

        self.fig = IntensityMapWidget(self.spectra, mode="powder")

        super().__init__(
            children=[
                self.settings_intensity,
                self.fig,
            ],
        )

    def _on_plot_button_clicked(self, change=None):
        self.parameters.update(
            {
                "weighting": self.children[0].weight_button.value,
                "q_spacing": self.children[0].slider_q_spacing.value,
                "energy_broadening": self.children[0].slider_energy_broadening.value,
                "ebins": self.children[0].slider_energy_bins.value,
                "temperature": self.children[0].slider_T.value,
                "q_min": self.children[0].slider_qmin.value,
                "q_max": self.children[0].slider_qmax.value,
                "npts": self.children[0].slider_npts.value,
            }
        )
        parameters = AttrDict(self.parameters)

        self.spectra, self.parameters = produce_powder_data(
            params=parameters, fc=self.fc, plot=False
        )
        self.fig._update_spectra(self.spectra, mode="powder")

    def download_data(self, _=None):
        filename = "powder_spectra.json"
        my_dict = self.spectra.to_dict()
        my_dict.update(
            {
                "weighting": self.children[0].weight_button.value,
                "q_spacing": self.children[0].slider_q_spacing.value,
                "energy_broadening": self.children[0].slider_energy_broadening.value,
                "ebins": self.children[0].slider_energy_bins.value,
                "temperature": self.children[0].slider_T.value,
                "q_min": self.children[0].slider_qmin.value,
                "q_max": self.children[0].slider_qmax.value,
                "npts": self.children[0].slider_npts.value,
            }
        )
        json_str = json.dumps(jsanitize(my_dict))
        b64_str = base64.b64encode(json_str.encode()).decode()
        self._download(payload=b64_str, filename=filename)

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


#### SETTINGS WIDGET:


class IntensitySettingsWidget(ipw.VBox):

    title_intensity = ipw.HTML("<h3>Dynamic structure factor</h3>")

    def __init__(self, mode="intensity", **kwargs):

        self.mode = mode

        self.title_intensity = ipw.HTML("<h3>Dynamic structure factor</h3>")

        self.slider_q_spacing = ipw.FloatSlider(
            min=0.01,
            max=1,
            step=0.01,
            value=0.01,
            description="&Delta;q (1/A)",
            tooltip="q spacing in 1/A",
        )

        self.slider_energy_broadening = ipw.FloatSlider(
            min=0.01,
            max=10,
            step=0.01,
            value=1,
            description="&Delta;E (meV)",
            tooltip="Energy broadening in meV",
        )

        self.slider_energy_bins = ipw.IntSlider(
            min=1,
            max=5000,
            step=1,
            value=1000,
            description="energy bins",
        )

        self.slider_T = ipw.IntSlider(
            min=0,
            max=1000,
            step=1,
            value=0,
            description="T (K)",
            disabled=False,
        )

        self.weight_button = ipw.ToggleButtons(
            options=[
                ("Coherent", "coherent"),
                ("DOS", "dos"),
            ],
            value="coherent",
            description="Weighting:",
            disabled=False,
            style={"description_width": "initial"},
        )

        self.plot_button = ipw.Button(
            description="Plot",
            icon="pencil",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(width="auto"),
        )

        self.reset_button = ipw.Button(
            description="Reset",
            icon="recycle",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(width="auto"),
        )

        self.download_button = ipw.Button(
            description="Download Data",
            icon="download",
            button_style="primary",
            disabled=True,
            layout=ipw.Layout(width="auto"),
        )

        self.reset_button.on_click(self._reset_settings)
        self.weight_button.observe(self._on_weight_button_change, names="value")

        if self.mode == "intensity":
            super().__init__(
                children=[
                    self.title_intensity,
                    self.slider_q_spacing,
                    self.slider_energy_broadening,
                    self.slider_energy_bins,
                    self.slider_T,
                    self.weight_button,
                    ipw.HBox(
                        [self.reset_button, self.plot_button, self.download_button]
                    ),
                ],
            )
        else:  # powder spectra. adding npts, q_min, q_max.
            self.title_intensity = ipw.HTML("<h3>Powder map</h3>")
            self.slider_qmin = ipw.FloatSlider(
                min=0,
                max=10,
                step=0.01,
                value=0,
                description="q<sub>min</sub> (1/A)",
            )

            self.slider_qmax = ipw.FloatSlider(
                min=0,
                max=10,
                step=0.01,
                value=1,
                description="q<sub>max</sub> (1/A)",
            )

            self.slider_npts = ipw.IntSlider(
                min=1,
                max=500,
                step=1,
                value=100,
                description="npts",
                tooltip="Number of points to be used in the average sphere.",
            )

            super().__init__(
                children=[
                    self.title_intensity,
                    self.slider_q_spacing,
                    self.slider_qmin,
                    self.slider_qmax,
                    self.slider_npts,
                    self.slider_energy_broadening,
                    self.slider_energy_bins,
                    self.slider_T,
                    self.weight_button,
                    ipw.HBox(
                        [self.reset_button, self.plot_button, self.download_button]
                    ),
                ],
            )

    def _reset_settings(self, _):
        self.slider_q_spacing.value = 0.01
        self.slider_energy_broadening.value = 1
        self.slider_energy_bins.value = 1000
        self.slider_T.value = 0
        self.weight_button.value = "coherent"

        if self.mode != "intensity":  # powder.
            self.slider_qmin.value = 0
            self.slider_qmax.value = 1
            self.slider_npts.value = 100

    def _on_weight_button_change(self, change):
        if change["new"] != change["old"]:
            self.slider_T.value = 0
            self.slider_T.disabled = True if change["new"] == "dos" else False
