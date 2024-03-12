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


def generate_force_constant_instance(
    phonopy_calc=None,
    path: str = None,
    summary_name: str = None,
    born_name: Optional[str] = None,
    fc_name: str = "FORCE_CONSTANTS",
    fc_format: Optional[str] = None,
):
    """
    Basically allows to obtain the ForceConstants instance from phonopy, both via files (from the second
    input parameters we have the same one of `euphonic.ForceConstants.from_phonopy`), or via a
    PhonopyCalculation instance. Respectively, the two ways will support independent euphonic app and integration
    of Euphonic into aiidalab.
    """
    blockPrint()

    ####### This is done to support the detached app (from aiidalab) with the same code:
    if path and summary_name:
        fc = euphonic.ForceConstants.from_phonopy(
            path=path,
            summary_name=summary_name,
            fc_name=fc_name,
        )
        return fc
    elif not phonopy_calc:
        raise NotImplementedError(
            "Please provide or the files or the phonopy calculation node."
        )

    ####### This is almost copied from PhonopyCalculation and is done to support functionalities in aiidalab env:
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
        # which is the case

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
        # each spectra has the .x_tick_labels attribute, for the bands.
        shift = False
        for k in j.x_tick_labels:
            ticks_positions.append(k[0])
            # ticks_labels.append("Gamma") if k[1] == '$\\Gamma$' else ticks_labels.append(k[1])
            ticks_labels.append(k[1])

            # Here below we check if we are starting a new group,
            # i.e. if the xticks count is starting again from 0
            # I also need to shift correctly the next index, which
            # refers to the zero of the ticks_positions[-1].
            if len(ticks_positions) > 1:
                if ticks_positions[-1] < ticks_positions[-2] or shift:
                    if ticks_positions[-1] == 0:  # new linear path

                        ticks_positions.pop()
                        last = ticks_labels.pop().strip()

                        # if the same index, do not join, just write once
                        if ticks_labels[-1].strip() != last:
                            ticks_labels[-1] = ticks_labels[-1].strip() + "|" + last
                        # the shift is needed because if this index was zero,
                        # the next one has to be shifted because it means that
                        # the index counting was restarted from zero,
                        # i.e. this is a new linear path.

                        shift = True
                    else:
                        ticks_positions[-1] = ticks_positions[-1] + ticks_positions[-2]

                if ticks_labels[-1] == ticks_labels[-2]:
                    ticks_positions.pop()
                    ticks_labels.pop()

    return final_xspectra, final_zspectra, ticks_positions, ticks_labels


# # Intensity map widget
class StructureFactorPlotWidget(ipw.VBox):
    """
    Widget to plot the Structure Factor for single crystals or powder samples.
    It takes as input the spectra as generated via the
    `produce_bands_weigthed_data` or `produce_powder_data` functions, called in the
    __init__ of the master widgets, respectively: `SingleCrystalFullWidget` and `PowderFullWidget`.

    NB: Intensity is relative to the maximum intensity at T=0K.
    """

    def __init__(self, spectra, mode="intensity", intensity_ref_0K=1, **kwargs):

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

        self.message_fig = ipw.HTML("")
        self.message_fig.layout.display = "none"

        self.fig = go.FigureWidget()

        self.fig.add_trace(
            go.Heatmap(
                z=final_zspectra.T / intensity_ref_0K,
                y=y[:, 0],
                x=None if mode == "intensity" else x[0],
                showscale=False,
            )
        )

        if self.fig.layout.images:
            for image in self.fig.layout.images:
                image["scl"] = 2  # Set the scale for each image

        if mode == "intensity":
            self.fig.update_layout(
                xaxis=dict(
                    tickmode="array", tickvals=ticks_positions, ticktext=ticks_labels
                )
            )
        else:
            self.fig["layout"]["xaxis"].update(title="|q| (1/A)")

        self.fig["layout"]["xaxis"].update(
            range=[min(final_xspectra), max(final_xspectra)]
        )
        self.fig["layout"]["yaxis"].update(range=[min(y[:, 0]), max(y[:, 0])])
        self.fig["layout"]["yaxis"].update(title="THz")

        self.fig.update_layout(
            height=450,
            width=650,
            margin=dict(l=20, r=20, t=20, b=20),
        )
        # Update x-axis and y-axis to enable autoscaling
        self.fig.update_xaxes(autorange=True)
        self.fig.update_yaxes(autorange=True)

        self.slider_intensity = ipw.FloatRangeSlider(
            value=[0, 100],  # Default selected interval
            min=0,
            max=100,
            step=0.1,
            description="Intensity window (%):",
            orientation="vertical",
            readout=True,
            readout_format=".0f",
            layout=ipw.Layout(
                width="20%",
            ),
        )

        self.slider_intensity.observe(self._update_intensity_filter, "value")

        # Create and show figure
        super().__init__(
            children=[self.message_fig, ipw.HBox([self.fig, self.slider_intensity])],
        )

    def _update_spectra(
        self,
        spectra,
        mode="intensity",
        intensity_ref_0K=1,
    ):

        if mode == "intensity":
            (
                final_xspectra,
                final_zspectra,
                ticks_positions,
                ticks_labels,
            ) = generated_curated_data(
                spectra,
            )
            # Data to contour is the sum of two Gaussian functions.
            x, y = np.meshgrid(spectra[0].x_data.magnitude, spectra[0].y_data.magnitude)
        else:
            final_zspectra = spectra.z_data.magnitude
            final_xspectra = spectra.x_data.magnitude
            # Data to contour is the sum of two Gaussian functions.
            x, y = np.meshgrid(spectra.x_data.magnitude, spectra.y_data.magnitude)

        # If I do this
        #   self.data = ()
        # I have a delay in the plotting, I have blank plot while it
        # is adding the new trace (see below); So, I will instead do the
        # re-assignement of the self.data = [self.data[1]] afterwards.

        x = None if mode == "intensity" else x[0]
        self.fig.add_trace(
            go.Heatmap(
                z=final_zspectra.T / intensity_ref_0K,
                y=y[:, 0],
                x=x,
                showscale=False,
            )
        )

        # change the path wants also a change in the labels.
        # this is delays things
        if mode == "intensity":
            self.fig.update_layout(
                xaxis=dict(
                    tickmode="array", tickvals=ticks_positions, ticktext=ticks_labels
                )
            )
        else:
            self.fig["layout"]["xaxis"].update(title="|q| (1/A)")

        self.fig.data = [self.fig.data[1]]

        # We should do a check, if we have few points (<200?) provide like a warning..
        # Also decise less than what, 30%, 50%...?
        visible_points = len(np.where(final_zspectra.T / intensity_ref_0K > 0.5)[0])
        if visible_points < 1000:
            message = f"Only {visible_points} points have intensity higher than 50%"
            self.message_fig.value = message
            self.message_fig.layout.display = "block"
        else:
            self.message_fig.layout.display = "none"

    def _update_intensity_filter(self, change):
        # the value of the intensity slider is in fractions of the max.
        if change["new"] != change["old"]:
            self.fig.data[0].zmax = (
                change["new"][1] * np.max(self.fig.data[0].z) / 100
            )  # above this, it is all yellow, i.e. max intensity.
            self.fig.data[0].zmin = (
                change["new"][0] * np.max(self.fig.data[0].z) / 100
            )  # below this, it is all blue, i.e. zero intensity


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

        self.settings_intensity = StructureFactorSettingsWidget()
        self.settings_intensity.plot_button.on_click(self._on_plot_button_clicked)
        self.settings_intensity.download_button.on_click(self.download_data)

        # This is used in order to have an overall intensity scale.
        self.intensity_ref_0K = np.max(self.spectra[0].z_data.magnitude)  # CHANGED

        self.map_widget = StructureFactorPlotWidget(
            self.spectra, intensity_ref_0K=self.intensity_ref_0K
        )  # CHANGED

        super().__init__(
            children=[
                self.settings_intensity.title_intensity,
                self.map_widget,
                self.settings_intensity,
            ],
        )

    def _on_plot_button_clicked(self, change=None):
        self.parameters.update(
            {
                "weighting": self.settings_intensity.weight_button.value,
                "q_spacing": self.settings_intensity.slider_q_spacing.value,
                "energy_broadening": self.settings_intensity.slider_energy_broadening.value,
                "ebins": self.settings_intensity.slider_energy_bins.value,
                "temperature": self.settings_intensity.slider_T.value,
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
        self.map_widget._update_spectra(
            self.spectra, intensity_ref_0K=self.intensity_ref_0K
        )  # CHANGED

    def download_data(self, _=None):
        """
        Download both the ForceConstants and the spectra json files.
        """
        force_constants_dict = fc.to_dict()

        filename = "single_crystal.json"
        my_dict = {}
        for branch in range(len(self.spectra)):
            my_dict[str(branch)] = self.spectra[branch].to_dict()
        my_dict.update(
            {
                "weighting": self.settings_intensity.weight_button.value,
                "q_spacing": self.settings_intensity.slider_q_spacing.value,
                "energy_broadening": self.settings_intensity.slider_energy_broadening.value,
                "ebins": self.settings_intensity.slider_energy_bins.value,
                "temperature": self.settings_intensity.slider_T.value,
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
        coordinates = []
        labels = []
        path = self.settings_intensity.custom_kpath_text.value
        linear_paths = path.split("|")
        for i in linear_paths:
            scoords = []
            s = i.split("-")
            for k in s:
                labels.append(k.strip())
                # AAA missing support for fractions.
                l = tuple([float(kk) for kk in k.strip().split(" ")])
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


class PowderFullWidget(ipw.VBox):
    """
    The Widget to display specifically the Intensity map of Dynamical structure
    factor for powder samples.

    The scattering lengths used in the `produce_bands_weigthed_data` function
    are tabulated (Euphonic/euphonic/data/sears-1992.json)
    and are from Sears (1992) Neutron News 3(3) pp26--37.
    """

    def __init__(self, fc, **kwargs):

        self.fc = fc

        self.spectra, self.parameters = produce_powder_data(
            params=parameters_powder, fc=self.fc
        )

        self.settings_intensity = StructureFactorSettingsWidget(mode="powder")
        self.settings_intensity.plot_button.on_click(self._on_plot_button_clicked)
        self.settings_intensity.download_button.on_click(self.download_data)

        # This is used in order to have an overall intensity scale.
        self.intensity_ref_0K = np.max(self.spectra.z_data.magnitude)  # CHANGED

        self.map_widget = StructureFactorPlotWidget(
            self.spectra, mode="powder", intensity_ref_0K=self.intensity_ref_0K
        )  # CHANGED

        super().__init__(
            children=[
                self.settings_intensity.title_intensity,
                self.map_widget,
                self.settings_intensity,
            ],
        )

    def _on_plot_button_clicked(self, change=None):
        self.parameters.update(
            {
                "weighting": self.settings_intensity.weight_button.value,
                "q_spacing": self.settings_intensity.slider_q_spacing.value,
                "energy_broadening": self.settings_intensity.slider_energy_broadening.value,
                "ebins": self.settings_intensity.slider_energy_bins.value,
                "temperature": self.settings_intensity.slider_T.value,
                "q_min": self.settings_intensity.slider_qmin.value,
                "q_max": self.settings_intensity.slider_qmax.value,
                "npts": self.settings_intensity.slider_npts.value,
            }
        )
        parameters_powder = AttrDict(self.parameters)

        self.spectra, self.parameters = produce_powder_data(
            params=parameters_powder, fc=self.fc, plot=False
        )

        self.settings_intensity.plot_button.disabled = True
        self.map_widget._update_spectra(
            self.spectra, mode="powder", intensity_ref_0K=self.intensity_ref_0K
        )  # CHANGED

    def download_data(self, _=None):
        """
        Download both the ForceConstants and the spectra json files.
        """
        force_constants_dict = fc.to_dict()

        filename = "powder"
        my_dict = self.spectra.to_dict()
        my_dict.update(
            {
                "weighting": self.settings_intensity.weight_button.value,
                "q_spacing": self.settings_intensity.slider_q_spacing.value,
                "energy_broadening": self.settings_intensity.slider_energy_broadening.value,
                "ebins": self.settings_intensity.slider_energy_bins.value,
                "temperature": self.settings_intensity.slider_T.value,
                "q_min": self.settings_intensity.slider_qmin.value,
                "q_max": self.settings_intensity.slider_qmax.value,
                # "npts": self.settings_intensity.slider_npts.value,
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


#### SETTINGS WIDGET:


class StructureFactorSettingsWidget(ipw.VBox):

    """
    Collects all the button and widget used to define settings for Neutron dynamic structure factor,
    both single crystal or powder.
    """

    title_intensity = ipw.HTML("<h3>Neutron dynamic structure factor</h3>")

    def __init__(self, mode="intensity", **kwargs):

        self.mode = mode

        self.title_intensity = ipw.HTML(
            "<h3>Neutron dynamic structure factor - Single Crystal</h3>"
        )
        self.specification_intensity = ipw.HTML(
            "(Intensity is relative to the maximum intensity at T=0K)"
        )

        self.slider_q_spacing = ipw.FloatSlider(
            min=0.005,
            max=0.5,
            step=0.005,
            value=0.001,
            description="q step (1/A)",
            tooltip="q spacing in 1/A",
        )
        self.slider_q_spacing.observe(self._on_setting_changed, names="value")

        self.slider_energy_broadening = ipw.FloatSlider(
            min=0.01,
            max=10,
            step=0.01,
            value=1,
            description="&Delta;E (THz)",
            tooltip="Energy broadening in THz",
        )
        self.slider_energy_broadening.observe(self._on_setting_changed, names="value")

        self.slider_energy_bins = ipw.IntSlider(
            min=1,
            max=5000,
            step=1,
            value=1000,
            description="#E bins",
            tooltip="Number of energy bins",
        )
        self.slider_energy_bins.observe(self._on_setting_changed, names="value")

        self.slider_T = ipw.IntSlider(
            min=0,
            max=1000,
            step=1,
            value=0,
            description="T (K)",
            disabled=False,
        )
        self.slider_T.observe(self._on_setting_changed, names="value")

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
        self.weight_button.observe(self._on_weight_button_change, names="value")

        self.custom_kpath_description = ipw.HTML(
            """<b>Custom q-points path for the structure factor</b>: <br>
            we can provide it via a specific format: <br>
            (1) each linear path should be divided by '|';<br>
            (2) each path is composed of 'qxi qyi qzi - qxf qyf qzf' where qxi and qxf are, respectively,
            the start and end x-coordinate of the q direction, in reciprocal lattice units (rlu).<br>
            An example path is: '0 0 0 - 1 1 1 | 1 1 1 - 0.5 0.5 0.5'. <br>
            For now, we do not support fractions (i.e. we accept 0.5 but not 1/2).
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

        # Please note: if you change the order of the widgets below, it will
        # affect the usage of the children[0] below.
        if self.mode == "intensity":
            super().__init__(
                children=[
                    # self.title_intensity,
                    ipw.HBox(
                        [self.reset_button, self.plot_button, self.download_button]
                    ),
                    ipw.HBox(
                        [
                            ipw.VBox(
                                [
                                    self.specification_intensity,
                                    self.slider_q_spacing,
                                    self.slider_energy_broadening,
                                    self.slider_energy_bins,
                                    self.slider_T,
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
                                    width="50%",
                                ),
                            ),
                        ],  # end of HBox children
                    ),
                ],
            )
        else:  # powder spectra. adding npts, q_min, q_max.
            self.title_intensity = ipw.HTML(
                "<h3>Neutron dynamic structure factor - Powder sample</h3>"
            )
            self.slider_qmin = ipw.FloatSlider(
                min=0,
                max=10,
                step=0.01,
                value=0,
                description="q<sub>min</sub> (1/A)",
            )
            self.slider_qmin.observe(self._on_setting_changed, names="value")

            self.slider_qmax = ipw.FloatSlider(
                min=0,
                max=10,
                step=0.01,
                value=1,
                description="q<sub>max</sub> (1/A)",
            )
            self.slider_qmax.observe(self._on_setting_changed, names="value")

            self.slider_npts = ipw.IntSlider(
                min=1,
                max=500,
                step=1,
                value=100,
                description="npts",
                tooltip="Number of points to be used in the average sphere.",
            )
            self.slider_npts.observe(self._on_setting_changed, names="value")

            super().__init__(
                children=[
                    # self.title_intensity,
                    ipw.HBox(
                        [self.reset_button, self.plot_button, self.download_button]
                    ),
                    ipw.HBox(
                        [
                            ipw.VBox(
                                [
                                    self.slider_q_spacing,
                                    self.slider_qmin,
                                    self.slider_qmax,
                                    # self.slider_npts,
                                ],
                                layout=ipw.Layout(
                                    width="50%",
                                ),
                            ),
                            ipw.VBox(
                                [
                                    self.slider_energy_broadening,
                                    self.slider_energy_bins,
                                    self.slider_T,
                                    self.weight_button,
                                ],
                                layout=ipw.Layout(
                                    width="50%",
                                ),
                            ),
                        ],  # end of HBox children
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
        else:
            self.custom_kpath_text.value = ""

    def _on_plot_button_changed(self, change):
        if change["new"] != change["old"]:
            self.download_button.disabled = not change["new"]

    def _on_weight_button_change(self, change):
        if change["new"] != change["old"]:
            self.slider_T.value = 0
            self.slider_T.disabled = True if change["new"] == "dos" else False
            self.plot_button.disabled = False

    def _on_setting_changed(self, change):
        self.plot_button.disabled = False


###### START for detached app:

# Upload buttons
class UploadPhonopyDataWidget(ipw.VBox):
    def __init__(self, **kwargs):

        self.upload_phonopy_button = ipw.FileUpload(
            description="upload phonopy YAML file",
            multiple=False,
            layout={"width": "initial"},
        )

        super().__init__(children=[self.upload_phonopy_button], **kwargs)

    def _read_phonopy_files(self, fname, content):
        suffix = "".join(pathlib.Path(fname).suffixes)

        with tempfile.NamedTemporaryFile(suffix=suffix) as temp_file:
            temp_file.write(content)
            temp_file.flush()
            try:
                if suffix == ".yaml":
                    fc = generate_force_constant_instance(
                        path=pathlib.Path(fname),
                        summary_name=temp_file.name,
                    )
                # else:
                #    structures = get_ase_from_file(temp_file.name)
            except ValueError as e:
                self._status_message.message = f"""
                    <div class="alert alert-danger">ERROR: {e}</div>
                    """
                return None
            except KeyError:
                self._status_message.message = f"""
                    <div class="alert alert-danger">ERROR: Could not parse file {fname}</div>
                    """
                return None

            return fc


#### END for detached app


##### START OVERALL WIDGET TO DISPLAY EVERYTHING:


class EuphonicSuperWidget(ipw.VBox):
    """
    Widget that will include everything,
    from the upload widget to the tabs with single crystal and powder predictions.
    """

    def __init__(self, mode="aiidalab-qe app plugin", fc=None):

        self.upload_widget = UploadPhonopyDataWidget()

        self.tab_widget = ipw.Tab()
        self.tab_widget.layout.display = "none"
        self.tab_widget.set_title(0, "Single crystal")
        self.tab_widget.set_title(1, "Powder sample")

        if mode == "aiidalab-qe app plugin":
            self.upload_widget.layout.display = "none"

            self.tab_widget.children.append(SingleCrystalFullWidget(fc))
            self.tab_widget.children.append(PowderFullWidget(fc))

            self.tab_widget.layout.display = "block"
        else:
            self.upload_widget.children[0].observe(self._on_upload_files, "value")

        super().__init__(
            children=[
                self.upload_widget,
                self.tab_widget,
            ],
        )

    def _on_upload_files(self, change):

        if change["new"] != change["old"]:

            self.tab_widget.children = ()

            fc = self.upload_widget._read_phonopy_files(
                "phonopy.yaml",
                self.upload_widget.children[0].value["phonopy.yaml"]["content"],
            )

            self.tab_widget.children = (
                SingleCrystalFullWidget(fc),
                PowderFullWidget(fc),
            )

            self.tab_widget.layout.display = "block"
