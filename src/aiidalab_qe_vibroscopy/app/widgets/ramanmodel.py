from __future__ import annotations
from aiidalab_qe.common.mvc import Model
import traitlets as tl
from aiida.common.extendeddicts import AttributeDict
from ase.atoms import Atoms
from IPython.display import display
import numpy as np
from aiida_vibroscopy.utils.broadenings import multilorentz
import plotly.graph_objects as go
import base64
import json


class RamanModel(Model):
    vibro = tl.Instance(AttributeDict, allow_none=True)
    input_structure = tl.Instance(Atoms, allow_none=True)
    spectrum_type = tl.Unicode()

    plot_type_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ("Powder", "powder"),
            ("Single Crystal", "single_crystal"),
        ],
    )

    plot_type = tl.Unicode("powder")
    temperature = tl.Float(300)
    frequency_laser = tl.Float(532)
    pol_incoming = tl.Unicode("0 0 1")
    pol_outgoing = tl.Unicode("0 0 1")
    broadening = tl.Float(10.0)
    separate_polarizations = tl.Bool(False)

    frequencies = []
    intensities = []

    frequencies_depolarized = []
    intensities_depolarized = []

    # Active modes
    active_modes_options = tl.List(
        trait=tl.Tuple((tl.Unicode(), tl.Int())), default_value=[]
    )
    active_mode = tl.Int()
    amplitude = tl.Float(3.0)

    supercell_0 = tl.Int(1)
    supercell_1 = tl.Int(1)
    supercell_2 = tl.Int(1)

    use_nac_direction = tl.Bool(False)
    nac_direction = tl.Unicode("0 0 1")

    def fetch_data(self):
        """Fetch the Raman data from the VibroWorkChain"""
        self.raman_data = self.get_vibrational_data(self.vibro)
        self.raw_frequencies, self.eigenvectors, self.labels = (
            self.raman_data.run_active_modes(
                selection_rule=self.spectrum_type.lower(),
            )
        )
        self.rounded_frequencies = [
            round(frequency, 3) for frequency in self.raw_frequencies
        ]
        self.active_modes_options = self._get_active_modes_options()
        self.active_mode = 0

    def _get_active_modes_options(self):
        active_modes_options = [
            (f"{index + 1}: {value}", index)
            for index, value in enumerate(self.rounded_frequencies)
        ]

        return active_modes_options

    def update_data(self):
        """
        Update the plot data based on the selected spectrum type, plot type, and configuration.
        """
        if self.plot_type == "powder":
            self._update_powder_data()
        else:
            self._update_single_crystal_data()

    def _update_powder_data(self):
        """
        Update data for the powder plot, handling both Raman and IR spectra.
        """
        dir_nac_direction, _ = self._check_inputs_correct(self.nac_direction)
        if self.spectrum_type == "Raman":
            (
                polarized_intensities,
                depolarized_intensities,
                frequencies,
                _,
            ) = self.raman_data.run_powder_raman_intensities(
                frequencies=self.frequency_laser,
                temperature=self.temperature,
                nac_direction=dir_nac_direction if self.use_nac_direction else None,
            )

            if self.separate_polarizations:
                self.frequencies, self.intensities = self.generate_plot_data(
                    frequencies,
                    polarized_intensities,
                    self.broadening,
                )
                self.frequencies_depolarized, self.intensities_depolarized = (
                    self.generate_plot_data(
                        frequencies,
                        depolarized_intensities,
                        self.broadening,
                    )
                )
            else:
                combined_intensities = polarized_intensities + depolarized_intensities
                self.frequencies, self.intensities = self.generate_plot_data(
                    frequencies,
                    combined_intensities,
                    self.broadening,
                )
                self.frequencies_depolarized, self.intensities_depolarized = [], []

        elif self.spectrum_type == "IR":
            (
                intensities,
                frequencies,
                _,
            ) = self.raman_data.run_powder_ir_intensities(
                nac_direction=dir_nac_direction if self.use_nac_direction else None,
            )
            self.frequencies, self.intensities = self.generate_plot_data(
                frequencies,
                intensities,
                self.broadening,
            )
            self.frequencies_depolarized, self.intensities_depolarized = [], []

    def _update_single_crystal_data(self):
        """
        Update data for the single crystal plot, handling both Raman and IR spectra.
        """
        dir_incoming, _ = self._check_inputs_correct(self.pol_incoming)
        dir_nac_direction, _ = self._check_inputs_correct(self.nac_direction)

        if self.spectrum_type == "Raman":
            dir_outgoing, _ = self._check_inputs_correct(self.pol_outgoing)
            (
                intensities,
                frequencies,
                _,
            ) = self.raman_data.run_single_crystal_raman_intensities(
                pol_incoming=dir_incoming,
                pol_outgoing=dir_outgoing,
                frequencies=self.frequency_laser,
                temperature=self.temperature,
                nac_direction=dir_nac_direction if self.use_nac_direction else None,
            )
        elif self.spectrum_type == "IR":
            (
                intensities,
                frequencies,
                _,
            ) = self.raman_data.run_single_crystal_ir_intensities(
                pol_incoming=dir_incoming,
                nac_direction=dir_nac_direction if self.use_nac_direction else None,
            )

        self.frequencies, self.intensities = self.generate_plot_data(
            frequencies, intensities
        )
        self.frequencies_depolarized, self.intensities_depolarized = [], []

    def update_plot(self, plot):
        """
        Update the Raman plot based on the selected plot type and configuration.

        Parameters:
            plot: The plotly.graph_objs.Figure widget to update.
        """
        update_function = (
            self._update_powder_plot
            if self.plot_type == "powder"
            else self._update_single_crystal_plot
        )
        update_function(plot)

    def _update_powder_plot(self, plot):
        """
        Update the powder Raman plot.

        Parameters:
            plot: The plotly.graph_objs.Figure widget to update.
        """
        if self.separate_polarizations:
            self._update_polarized_and_depolarized(plot)
        else:
            self._clear_depolarized_and_update(plot)

    def _update_polarized_and_depolarized(self, plot):
        """
        Update the plot when polarized and depolarized data are separate.

        Parameters:
            plot: The plotly.graph_objs.Figure widget to update.
        """
        if len(plot.data) == 1:
            self._update_trace(
                plot.data[0], self.frequencies, self.intensities, "Polarized"
            )
            plot.add_trace(
                go.Scatter(
                    x=self.frequencies_depolarized,
                    y=self.intensities_depolarized,
                    name="Depolarized",
                )
            )
            plot.layout.title.text = f"Powder {self.spectrum_type} spectrum"
        elif len(plot.data) == 2:
            self._update_trace(
                plot.data[0], self.frequencies, self.intensities, "Polarized"
            )
            self._update_trace(
                plot.data[1],
                self.frequencies_depolarized,
                self.intensities_depolarized,
                "Depolarized",
            )
            plot.layout.title.text = f"Powder {self.spectrum_type} spectrum"

    def _clear_depolarized_and_update(self, plot):
        """
        Clear depolarized data and update the plot.

        Parameters:
            plot: The plotly.graph_objs.Figure widget to update.
        """
        if len(plot.data) == 2:
            self._update_trace(plot.data[0], self.frequencies, self.intensities, "")
            plot.data[1].x = []
            plot.data[1].y = []
            plot.layout.title.text = f"Powder {self.spectrum_type} spectrum"
        elif len(plot.data) == 1:
            self._update_trace(plot.data[0], self.frequencies, self.intensities, "")
            plot.layout.title.text = f"Powder {self.spectrum_type} spectrum"

    def _update_single_crystal_plot(self, plot):
        """
        Update the single crystal Raman plot.

        Parameters:
            plot: The plotly.graph_objs.Figure widget to update.
        """
        if len(plot.data) == 2:
            self._update_trace(plot.data[0], self.frequencies, self.intensities, "")
            plot.data[1].x = []
            plot.data[1].y = []
            plot.layout.title.text = f"Single crystal {self.spectrum_type} spectrum"
        elif len(plot.data) == 1:
            self._update_trace(plot.data[0], self.frequencies, self.intensities, "")
            plot.layout.title.text = f"Single crystal {self.spectrum_type} spectrum"

    def _update_trace(self, trace, x_data, y_data, name):
        """
        Helper function to update a single trace in the plot.

        Parameters:
            trace: The trace to update.
            x_data: The new x-axis data.
            y_data: The new y-axis data.
            name: The name of the trace.
        """
        trace.x = x_data
        trace.y = y_data
        trace.name = name

    def get_vibrational_data(self, node):
        """
        Extract vibrational data from an IRamanWorkChain or HarmonicWorkChain node.

        Parameters:
            node: The workchain node containing IRaman or Harmonic data.

        Returns:
            The vibrational accuracy data (vibro) or None if not available.
        """
        # Determine the output node
        output_node = getattr(node, "iraman", None) or getattr(node, "harmonic", None)
        if not output_node:
            return None

        # Check for vibrational data and extract accuracy
        vibrational_data = getattr(output_node, "vibrational_data", None)
        if not vibrational_data:
            return None

        # Extract vibrational accuracy (prefer numerical_accuracy_4 if available)
        vibro = getattr(vibrational_data, "numerical_accuracy_4", None) or getattr(
            vibrational_data, "numerical_accuracy_2", None
        )

        return vibro

    def _check_inputs_correct(self, polarization):
        # Check if the polarization vectors are correct
        input_text = polarization
        input_values = input_text.split()
        dir_values = []
        if len(input_values) == 3:
            try:
                dir_values = [float(i) for i in input_values]
                return dir_values, True
            except:  # noqa: E722
                return dir_values, False
        else:
            return dir_values, False

    def generate_plot_data(
        self,
        frequencies: list[float],
        intensities: list[float],
        broadening: float = 10.0,
        x_range: list[float] | str = "auto",
        broadening_function=multilorentz,
        normalize: bool = True,
    ):
        frequencies = np.array(frequencies)
        intensities = np.array(intensities)

        if x_range == "auto":
            xi = max(0, frequencies.min() - 200)
            xf = frequencies.max() + 200
            x_range = np.arange(xi, xf, 1.0)

        y_range = broadening_function(x_range, frequencies, intensities, broadening)

        if normalize:
            y_range /= y_range.max()

        return x_range, y_range

    def modes_table(self):
        """Display table with the active modes."""
        # Create an HTML table with the active modes
        table_data = [list(x) for x in zip(self.rounded_frequencies, self.labels)]
        table_html = "<table>"
        table_html += "<tr><th>Frequencies (cm<sup>-1</sup>) </th><th> Label</th></tr>"
        for row in table_data:
            table_html += "<tr>"
            for cell in row:
                table_html += "<td style='text-align:center;'>{}</td>".format(cell)
            table_html += "</tr>"
        table_html += "</table>"

        return table_html

    def set_vibrational_mode_animation(self, weas):
        eigenvector = self.eigenvectors[self.active_mode]
        phonon_setting = {
            "eigenvectors": np.array(
                [[[real_part, 0] for real_part in row] for row in eigenvector]
            ),
            "kpoint": [0, 0, 0],  # optional
            "amplitude": self.amplitude,
            "factor": self.amplitude * 0.6,
            "nframes": 20,
            "repeat": [
                self.supercell_0,
                self.supercell_1,
                self.supercell_2,
            ],
            "color": "black",
            "radius": 0.1,
        }
        weas.avr.phonon_setting = phonon_setting
        return weas

    def download_data(self, _=None):
        filename = "spectra.json"
        if self.separate_polarizations:
            my_dict = {
                "Frequencies cm-1": self.frequencies.tolist(),
                "Polarized intensities": self.intensities.tolist(),
                "Depolarized intensities": self.intensities_depolarized.tolist(),
            }
        else:
            my_dict = {
                "Frequencies cm-1": self.frequencies.tolist(),
                "Intensities": self.intensities.tolist(),
            }
        json_str = json.dumps(my_dict)
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
            """.format(payload=payload, filename=filename)
        )
        display(javas)
