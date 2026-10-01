from __future__ import annotations
from aiidalab_qe.common.mvc import Model
import traitlets as tl
from aiida import orm
from ase.atoms import Atoms
from IPython.display import display
import numpy as np
from aiida_vibroscopy.utils.broadenings import multilorentz
import plotly.graph_objects as go
import base64
import json
from scipy.integrate import dblquad
from aiida_vibroscopy.utils.spectra import raman_prefactor
from aiidalab_qe_vibroscopy.utils.atom_selection import parse_atom_selection
from aiidalab_qe_vibroscopy.utils.mode_projection import (
    atom_participation,
    unitcell_to_primitive,
)


class RamanModel(Model):
    vibrational_data_uuid = tl.Unicode()
    selected_atoms = tl.Unicode()
    selected_only = tl.Bool(False)
    input_structure = tl.Instance(Atoms, allow_none=True)
    spectrum_type = tl.Unicode()

    plot_type_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ("Powder", "powder"),
            ("Single Crystal", "single_crystal"),
            ("2D average", "plane_average"),
        ],
    )

    plot_type = tl.Unicode("powder")

    plane_type_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ("XY", "xy"),
            ("YZ", "yz"),
            ("XZ", "xz"),
        ],
    )
    plane_type = tl.Unicode("xy")

    temperature = tl.Float(300)
    frequency_laser = tl.Float(532)
    pol_incoming = tl.Unicode("0 0 1")
    pol_outgoing = tl.Unicode("0 0 1")
    broadening = tl.Float(10.0)
    separate_polarizations = tl.Bool(False)

    frequencies = []
    intensities = []

    raw_frequencies = []
    raw_intensities = []
    raw_pol_intensities = []
    raw_depol_intensities = []

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

    def get_raman_data(self):
        """Load only for the immediate operation; retain the UUID in the model."""
        return orm.load_node(self.vibrational_data_uuid)

    def fetch_data(self):
        self._refresh_modes()
        self.active_mode = 0

    def _refresh_modes(self):
        direction, _ = self._check_inputs_correct(self.nac_direction)
        direction = direction if self.use_nac_direction else None
        key = (self.vibrational_data_uuid, self.spectrum_type, tuple(direction or []))
        if key == getattr(self, "_mode_key", None):
            return
        data = self.get_raman_data()
        frequencies, displacements, labels = data.run_active_modes(
            selection_rule=self.spectrum_type.lower(),
            nac_direction=direction,
        )
        if not len(frequencies):
            raise ValueError("No active modes are available for this NAC direction.")
        phonopy = data.get_phonopy_instance()
        mapping = unitcell_to_primitive(phonopy, self.input_structure)
        self.eigenvectors = np.asarray(displacements)[:, mapping, :]
        self._atom_participation = atom_participation(
            self.eigenvectors, self.input_structure.get_masses(), frequencies
        )
        self._mode_frequencies = np.asarray(frequencies)
        self.raw_frequencies = np.asarray(frequencies)
        self.labels = list(labels)
        self.rounded_frequencies = np.round(frequencies, 3).tolist()
        self.active_modes_options = self._get_active_modes_options()
        self.active_mode = min(self.active_mode, len(frequencies) - 1)
        self._mode_key = key

    def _get_active_modes_options(self):
        active_modes_options = [
            (f"{index + 1}: {value}", index)
            for index, value in enumerate(self.rounded_frequencies)
        ]

        return active_modes_options

    def _update_spectrum_options(self):
        if self.spectrum_type == "Raman":
            self.plot_type_options = [
                ("Powder", "powder"),
                ("Single Crystal", "single_crystal"),
                ("2D average", "plane_average"),
            ]
        else:
            self.plot_type_options = [
                ("Powder", "powder"),
                ("Single Crystal", "single_crystal"),
            ]

    def update_data(self):
        """
        Update the plot data based on the selected spectrum type, plot type, and configuration.
        """
        self.selected_indices = parse_atom_selection(
            self.selected_atoms,
            len(self.input_structure) if self.input_structure is not None else 0,
        )
        if self.vibrational_data_uuid:
            self._refresh_modes()
        if self.plot_type == "powder":
            self._update_powder_data()
        elif self.plot_type == "single_crystal":
            self._update_single_crystal_data()
        elif self.plot_type == "plane_average":
            self._update_plane_average_data()
        self._update_atom_projection()

    def _update_atom_projection(self):
        self.projected_intensities = np.array([])
        self.projected_depolarized = np.array([])
        self.projection_weights = np.array([])
        if not self.selected_indices:
            return
        if len(self._mode_frequencies) != len(self.raw_frequencies) or not np.allclose(
            self._mode_frequencies, self.raw_frequencies, atol=1e-5, rtol=1e-8
        ):
            raise ValueError("The participation modes do not match this spectrum.")
        self.projection_weights = self._atom_participation[
            :, self.selected_indices
        ].sum(axis=1)
        if self._has_separate_polarizations():
            self.projected_intensities = self._project_spectrum(
                self.raw_pol_intensities
            )
            self.projected_depolarized = self._project_spectrum(
                self.raw_depol_intensities
            )
        else:
            self.projected_intensities = self._project_spectrum(self.raw_intensities)

    def _has_separate_polarizations(self):
        return (
            self.spectrum_type == "Raman"
            and self.plot_type == "powder"
            and self.separate_polarizations
        )

    def _project_spectrum(self, intensities):
        _, total = self.generate_plot_data(
            self.raw_frequencies,
            intensities,
            self.broadening,
            x_range=self.frequencies,
            normalize=False,
        )
        _, selected = self.generate_plot_data(
            self.raw_frequencies,
            np.asarray(intensities) * self.projection_weights,
            self.broadening,
            x_range=self.frequencies,
            normalize=False,
        )
        scale = total.max(initial=0)
        return selected / scale if scale > 0 else np.zeros_like(selected)

    def _update_powder_data(self):
        """
        Update data for the powder plot, handling both Raman and IR spectra.
        """
        dir_nac_direction, _ = self._check_inputs_correct(self.nac_direction)
        if self.spectrum_type == "Raman":
            (
                self.raw_pol_intensities,
                self.raw_depol_intensities,
                self.raw_frequencies,
                _,
            ) = self.get_raman_data().run_powder_raman_intensities(
                frequency_laser=self.frequency_laser,
                temperature=self.temperature,
                nac_direction=dir_nac_direction if self.use_nac_direction else None,
            )

            if self.separate_polarizations:
                self.frequencies, self.intensities = self.generate_plot_data(
                    self.raw_frequencies,
                    self.raw_pol_intensities,
                    self.broadening,
                )
                self.frequencies_depolarized, self.intensities_depolarized = (
                    self.generate_plot_data(
                        self.raw_frequencies,
                        self.raw_depol_intensities,
                        self.broadening,
                    )
                )
            else:
                self.raw_intensities = (
                    self.raw_pol_intensities + self.raw_depol_intensities
                )
                self.frequencies, self.intensities = self.generate_plot_data(
                    self.raw_frequencies,
                    self.raw_intensities,
                    self.broadening,
                )
                self.frequencies_depolarized, self.intensities_depolarized = [], []
                self.raw_pol_intensities, self.raw_depol_intensities = [], []

        elif self.spectrum_type == "IR":
            (
                self.raw_intensities,
                self.raw_frequencies,
                _,
            ) = self.get_raman_data().run_powder_ir_intensities(
                nac_direction=dir_nac_direction if self.use_nac_direction else None,
            )
            self.frequencies, self.intensities = self.generate_plot_data(
                self.raw_frequencies,
                self.raw_intensities,
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
                self.raw_intensities,
                self.raw_frequencies,
                _,
            ) = self.get_raman_data().run_single_crystal_raman_intensities(
                pol_incoming=dir_incoming,
                pol_outgoing=dir_outgoing,
                frequency_laser=self.frequency_laser,
                temperature=self.temperature,
                nac_direction=dir_nac_direction if self.use_nac_direction else None,
            )
        elif self.spectrum_type == "IR":
            (
                self.raw_intensities,
                self.raw_frequencies,
                _,
            ) = self.get_raman_data().run_single_crystal_ir_intensities(
                pol_incoming=dir_incoming,
                nac_direction=dir_nac_direction if self.use_nac_direction else None,
            )

        self.frequencies, self.intensities = self.generate_plot_data(
            self.raw_frequencies, self.raw_intensities, self.broadening
        )
        self.frequencies_depolarized, self.intensities_depolarized = [], []

    def _update_plane_average_data(self):
        "Average the Raman susceptibility tensors over a plane."

        dir_nac_direction, _ = self._check_inputs_correct(self.nac_direction)

        def intensity(a, b, c, d):
            return dblquad(
                lambda t, x: np.abs(
                    a * np.cos(t) * np.cos(t + x)
                    + b * np.sin(t) * np.cos(t + x)
                    + c * np.cos(t) * np.sin(t + x)
                    + d * np.sin(t) * np.sin(t + x)
                )
                ** 2,
                0,
                2 * np.pi,
                lambda x: 0,
                lambda x: 2 * np.pi,
            )

        def get_plane_subtensor(raman_susc_tensor, plane):
            if plane == "xy":
                return raman_susc_tensor[np.ix_([0, 1], [0, 1])]
            elif plane == "yz":
                return raman_susc_tensor[np.ix_([1, 2], [1, 2])]
            elif plane == "xz":
                return raman_susc_tensor[np.ix_([0, 2], [0, 2])]

        raman_susc_tensor, self.raw_frequencies, _ = (
            self.get_raman_data().run_raman_susceptibility_tensors(
                nac_direction=dir_nac_direction if self.use_nac_direction else None,
            )
        )

        # Average the susceptibility tensors over frequencies at given plane
        intensities_plane = []

        for tensor in raman_susc_tensor:
            plane_subtensor = get_plane_subtensor(tensor, self.plane_type)
            a, b, c, d = plane_subtensor.flatten()
            avg_intensity, _ = intensity(a, b, c, d)
            intensities_plane.append(avg_intensity)

        intensities = np.array(intensities_plane)
        self.raw_intensities = intensities * raman_prefactor(
            frequency=self.raw_frequencies,
            frequency_laser=self.frequency_laser,
            temperature=self.temperature,
            absolute=True,
        )

        self.frequencies, self.intensities = self.generate_plot_data(
            self.raw_frequencies,
            self.raw_intensities,
            self.broadening,
        )

    def update_plot(self, plot):
        """Show the total and selected participation on the same intensity scale."""
        separate = self._has_separate_polarizations()
        channels = [
            (
                self.intensities,
                self.projected_intensities,
                "Polarized" if separate else "Total",
                "#555555",
                "#d62728",
            )
        ]
        if separate:
            channels.append(
                (
                    self.intensities_depolarized,
                    self.projected_depolarized,
                    "Depolarized",
                    "#999999",
                    "#ff7f0e",
                )
            )
        with plot.batch_update():
            plot.data = ()
            for total, selected, name, total_color, selected_color in channels:
                if not (self.selected_indices and self.selected_only):
                    plot.add_trace(
                        go.Scatter(
                            x=self.frequencies,
                            y=total,
                            name=name,
                            line={"width": 1.5, "color": total_color},
                        )
                    )
                if self.selected_indices:
                    label = "Selected-atom mode participation"
                    if separate:
                        label += f" ({name.lower()})"
                    plot.add_trace(
                        go.Scatter(
                            x=self.frequencies,
                            y=selected,
                            name=label,
                            line={"width": 3.5, "color": selected_color},
                        )
                    )
            title = {
                "powder": "Powder",
                "single_crystal": "Single crystal",
                "plane_average": f"{self.plane_type.upper()} plane average",
            }[self.plot_type]
            plot.layout.title.text = f"{title} {self.spectrum_type} spectrum"

    @staticmethod
    def get_vibrational_data(node):
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

        if isinstance(x_range, str) and x_range == "auto":
            xi = max(0, frequencies.min() - 200)
            xf = frequencies.max() + 200
            x_range = np.arange(xi, xf, 1.0)

        y_range = broadening_function(x_range, frequencies, intensities, broadening)

        if normalize:
            scale = y_range.max(initial=0)
            if scale > 0:
                y_range /= scale

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
        weas._widget.viewerStyle = {"width": "800px", "height": "600px"}
        weas.avr.phonon_setting = phonon_setting
        return weas

    def download_data(self, _=None):
        filename = "spectra.json"
        if self._has_separate_polarizations():
            my_dict = {
                "Frequencies cm-1": self.frequencies.tolist(),
                "Polarized intensities": self.intensities.tolist(),
                "Depolarized intensities": self.intensities_depolarized.tolist(),
                "Eigenvectors": self.eigenvectors.tolist(),
                "Raw Frequencies cm-1": self.raw_frequencies.tolist(),
                "Raw Intensities Polarized": self.raw_pol_intensities.tolist(),
                "Raw Intensities Depolarized": self.raw_depol_intensities.tolist(),
                "Labels": self.labels,
            }
        else:
            my_dict = {
                "Frequencies cm-1": self.frequencies.tolist(),
                "Intensities": self.intensities.tolist(),
                "Eigenvectors": self.eigenvectors.tolist(),
                "Raw Frequencies cm-1": self.raw_frequencies.tolist(),
                "Raw Intensities": self.raw_intensities.tolist(),
                "Labels": self.labels,
            }
        if self.selected_indices:
            my_dict.update(
                {
                    "Selected atoms (1-based)": [i + 1 for i in self.selected_indices],
                    "Projection definition": (
                        "Mass-weighted mode participation, averaged within degenerate "
                        "subspaces; shared normalization with the total spectrum."
                    ),
                    "Mode participation": self.projection_weights.tolist(),
                    "Selected-atom mode participation": self.projected_intensities.tolist(),
                }
            )
            if self._has_separate_polarizations():
                my_dict["Selected depolarized participation"] = (
                    self.projected_depolarized.tolist()
                )
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
