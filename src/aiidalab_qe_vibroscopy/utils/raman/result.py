"""Bands results view widgets

"""
from __future__ import annotations


from aiidalab_qe.common.panel import ResultPanel
import ipywidgets as ipw
import numpy as np
from IPython.display import HTML, clear_output, display
import base64
import json
from ase.visualize import view
import nglview as nv
from ase import Atoms


from aiida_vibroscopy.utils.broadenings import multilorentz


def plot_powder(
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


def export_iramanworkchain_data(node):

    """
    We have multiple choices: IR, RAMAN.
    """

    if not "vibronic" in node.outputs:
        return None
    else:
        if "iraman" in node.outputs.vibronic:
            output_node = node.outputs.vibronic.iraman
        elif "harmonic" in node.outputs.vibronic:
            output_node = node.outputs.vibronic.harmonic
        else:
            # we have raman and ir only if we run IRamanWorkChain or HarmonicWorkChain
            return None

    if "vibrational_data" in output_node:

        # We enable the possibility to provide both spectra.
        # We give as output or string, or the output node.

        spectra_data = {
            "Raman": None,
            "Ir": None,
        }

        vibrational_data = output_node.vibrational_data
        vibro = (
            vibrational_data.numerical_accuracy_4
            if hasattr(vibrational_data, "numerical_accuracy_4")
            else vibrational_data.numerical_accuracy_2
        )

        if "born_charges" in vibro.get_arraynames():
            (
                polarized_intensities,
                frequencies,
                labels,
            ) = vibro.run_powder_ir_intensities()
            total_intensities = polarized_intensities

            # sometimes IR/Raman has not active peaks by symmetry, or due to the fact that 1st order cannot capture them
            if len(total_intensities) == 0:
                spectra_data[
                    "Ir"
                ] = "No IR modes detected."  # explanation added in the main results script of the app.
            else:
                spectra_data["Ir"] = output_node

        if "raman_tensors" in vibro.get_arraynames():

            (
                polarized_intensities,
                depolarized_intensities,
                frequencies,
                labels,
            ) = vibro.run_powder_raman_intensities(frequency_laser=532, temperature=300)
            total_intensities = polarized_intensities + depolarized_intensities

            # sometimes IR/Raman has not active peaks by symmetry, or due to the fact that 1st order cannot capture them
            if len(total_intensities) == 0:
                spectra_data[
                    "Raman"
                ] = "No Raman modes detected."  # explanation added in the main results script of the app.
            else:
                spectra_data["Raman"] = output_node

        return spectra_data
    else:
        return None


class SpectrumPlotWidget(ipw.VBox):
    """Widget that allows different options for plotting Raman or IR Spectrum."""

    def __init__(self, node, output_node, spectrum_type, **kwargs):
        self.node = node
        self.output_node = output_node
        self.spectrum_type = spectrum_type

        # VibrationalData
        vibrational_data = self.output_node.vibrational_data
        self.vibro = (
            vibrational_data.numerical_accuracy_4
            if hasattr(vibrational_data, "numerical_accuracy_4")
            else vibrational_data.numerical_accuracy_2
        )

        self.description = ipw.HTML(
            f"""<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
            Select the type of {self.spectrum_type} spectrum to plot.
            </div>"""
        )

        self._plot_type = ipw.ToggleButtons(
            options=[
                ("Powder", "powder"),
                ("Single Crystal", "single_crystal"),
            ],
            value="powder",
            description="Spectrum type:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.temperature = ipw.FloatText(
            value=298.0,
            description="Temperature (K):",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.frequency_laser = ipw.FloatText(
            value=532.0,
            description="Laser frequency (nm):",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.pol_incoming = ipw.Text(
            value="0 0 1",
            description="Incoming polarization:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.pol_outgoing = ipw.Text(
            value="0 0 1",
            description="Outgoing polarization:",
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
        self.download_button = ipw.Button(
            description="Download Data",
            icon="download",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(width="auto", visibility="hidden"),
        )
        self.wrong_syntax = ipw.HTML(
            value="""<i class="fa fa-times" style="color:red;font-size:2em;" ></i> wrong syntax""",
            layout={"visibility": "hidden"},
        )
        self.broadening = ipw.FloatText(
            value=10.0,
            description="Broadening (cm-1):",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.spectrum_widget = ipw.Output()
        self.frequencies = []
        self.intensities = []
        self.polarization_out = ipw.Output()

        def download_data(_=None):
            filename = "spectra.json"
            my_dict = {
                "Frequencies cm-1": self.frequencies.tolist(),
                "Intensities": self.intensities.tolist(),
            }
            json_str = json.dumps(my_dict)
            b64_str = base64.b64encode(json_str.encode()).decode()
            self._download(payload=b64_str, filename=filename)

        # hide the outgoing pol, frequency and temperature if ir:
        if self.spectrum_type == "IR":
            self.temperature.layout.display = "none"
            self.frequency_laser.layout.display = "none"
            self.pol_outgoing.layout.display == "none"

        self._plot_type.observe(self._on_plot_type_change, names="value")
        self.plot_button.on_click(self._on_plot_button_clicked)
        self.download_button.on_click(download_data)
        super().__init__(
            children=[
                self.description,
                self._plot_type,
                self.temperature,
                self.frequency_laser,
                self.broadening,
                self.polarization_out,
                ipw.HBox([self.plot_button, self.download_button]),
                self.wrong_syntax,
                self.spectrum_widget,
            ]
        )

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

    def _on_plot_type_change(self, change):
        if change["new"] == "single_crystal":
            with self.polarization_out:
                clear_output()
                display(self.pol_incoming)
                if self.spectrum_type == "Raman":
                    display(self.pol_outgoing)
        else:
            self.pol_incoming.value = "0 0 1"
            self.pol_outgoing.value = "0 0 1"
            self.wrong_syntax.layout.visibility = "hidden"
            with self.polarization_out:
                clear_output()

    def _on_plot_button_clicked(self, change):
        if self._plot_type.value == "powder":
            # Powder spectrum
            if self.spectrum_type == "Raman":
                (
                    polarized_intensities,
                    depolarized_intensities,
                    frequencies,
                    labels,
                ) = self.vibro.run_powder_raman_intensities(
                    frequencies=self.frequency_laser.value,
                    temperature=self.temperature.value,
                )
                total_intensities = polarized_intensities + depolarized_intensities
            else:
                (
                    polarized_intensities,
                    frequencies,
                    labels,
                ) = self.vibro.run_powder_ir_intensities()
                total_intensities = polarized_intensities

            self.frequencies, self.intensities = plot_powder(
                frequencies, total_intensities, self.broadening.value,
            )
            self._display_figure()

        else:
            # Single crystal spectrum
            dir_incoming, correct_syntax_incoming = self._check_inputs_correct(
                self.pol_incoming
            )
            dir_outgoing, correct_syntax_outgoing = self._check_inputs_correct(
                self.pol_outgoing
            )
            if not correct_syntax_incoming or not correct_syntax_outgoing:
                self.wrong_syntax.layout.visibility = "visible"
                return
            else:
                self.wrong_syntax.layout.visibility = "hidden"
                if self.spectrum_type == "Raman":
                    (
                        intensities,
                        frequencies,
                        labels,
                    ) = self.vibro.run_single_crystal_raman_intensities(
                        pol_incoming=dir_incoming,
                        pol_outgoing=dir_incoming,
                        frequencies=self.frequency_laser.value,
                        temperature=self.temperature.value,
                    )
                else:
                    (
                        intensities,
                        frequencies,
                        labels,
                    ) = self.vibro.run_single_crystal_ir_intensities(
                        pol_incoming=dir_incoming
                    )

                self.frequencies, self.intensities = plot_powder(
                    frequencies, intensities
                )
                self._display_figure()

    def _check_inputs_correct(self, polarization):
        # Check if the polarization vectors are correct
        input_text = polarization.value
        input_values = input_text.split()
        dir_values = []
        if len(input_values) == 3:
            try:
                dir_values = [float(i) for i in input_values]
                return dir_values, True
            except:
                return dir_values, False
        else:
            return dir_values, False

    def _spectrum_widget(self):
        import plotly.graph_objects as go

        fig = go.FigureWidget(
            layout=go.Layout(
                title=dict(text=f"{self.spectrum_type} spectrum"),
                barmode="overlay",
            )
        )
        fig.layout.xaxis.title = "Wavenumber (cm-1)"
        fig.layout.yaxis.title = "Intensity (arb. units)"
        fig.layout.xaxis.nticks = 0
        fig.add_scatter(x=self.frequencies, y=self.intensities, name=f"")
        fig.update_layout(
            height=500,
            width=700,
            plot_bgcolor="white",
        )
        return fig

    def _display_figure(self):
        with self.spectrum_widget:
            clear_output()
            display(self._spectrum_widget())
            self.download_button.layout.visibility = "visible"


class ActiveModesWidget(ipw.VBox):
    """Widget that display an animation (nglview) of the active modes."""

    def __init__(self, node, output_node, spectrum_type, **kwargs):
        self.node = node
        self.output_node = output_node
        self.spectrum_type = spectrum_type

        # VibrationalData
        vibrational_data = self.output_node.vibrational_data
        self.vibro = (
            vibrational_data.numerical_accuracy_4
            if hasattr(vibrational_data, "numerical_accuracy_4")
            else vibrational_data.numerical_accuracy_2
        )

        # Raman or IR active modes
        selection_rule = self.spectrum_type.lower()
        frequencies, self.eigenvectors, self.labels = self.vibro.run_active_modes(
            selection_rule=selection_rule,
        )
        self.rounded_frequencies = [round(frequency, 3) for frequency in frequencies]

        # StructureData
        self.structure_ase = self.node.inputs.structure.get_ase()

        self.ngl = nv.NGLWidget()

        modes_values = [
            f"{index + 1}: {value}"
            for index, value in enumerate(self.rounded_frequencies)
        ]
        # Create Raman modes widget
        self.active_modes = ipw.Dropdown(
            options=modes_values,
            value=modes_values[0],  # Default value
            description="Select mode:",
            style={"description_width": "initial"},
        )

        self.amplitude = ipw.FloatText(
            value=3.0,
            description="Amplitude :",
            disabled=False,
            style={"description_width": "initial"},
        )

        self._supercell = [
            ipw.BoundedIntText(value=1, min=1, layout={"width": "40px"}),
            ipw.BoundedIntText(value=1, min=1, layout={"width": "40px"}),
            ipw.BoundedIntText(value=1, min=1, layout={"width": "40px"}),
        ]

        self.supercell_selector = ipw.HBox(
            [
                ipw.HTML(
                    description="Super cell:", style={"description_width": "initial"}
                )
            ]
            + self._supercell
        )

        self.modes_table = ipw.Output()
        self.animation = ipw.Output()
        self._display_table()
        self._select_active_mode(None)
        widget_list = [
            self.active_modes,
            self.amplitude,
            self._supercell[0],
            self._supercell[1],
            self._supercell[2],
        ]
        for elem in widget_list:
            elem.observe(self._select_active_mode, names="value")

        super().__init__(
            children=[
                ipw.HBox(
                    [
                        ipw.VBox(
                            [
                                ipw.HTML(
                                    value=f"<b> {self.spectrum_type} Active Modes </b>"
                                ),
                                self.modes_table,
                            ]
                        ),
                        ipw.VBox(
                            [
                                self.active_modes,
                                self.amplitude,
                                self.supercell_selector,
                                self.animation,
                            ],
                            layout={"justify_content": "center"},
                        ),
                    ]
                ),
            ]
        )

    def _display_table(self):
        """Display table with the active modes."""
        # Create an HTML table with the active modes
        table_data = [list(x) for x in zip(self.rounded_frequencies, self.labels)]
        table_html = "<table>"
        table_html += "<tr><th>Frequencies (cm-1) </th><th> Label</th></tr>"
        for row in table_data:
            table_html += "<tr>"
            for cell in row:
                table_html += "<td style='text-align:center;'>{}</td>".format(cell)
            table_html += "</tr>"
        table_html += "</table>"
        # Set layout to a fix size
        self.modes_table.layout = {
            "overflow": "auto",
            "height": "200px",
            "width": "150px",
        }
        with self.modes_table:
            clear_output()
            display(HTML(table_html))

    def _select_active_mode(self, change):
        """Display animation of the selected active mode."""
        self._animation_widget()
        with self.animation:
            clear_output()
            display(self.ngl)

    def _animation_widget(self):
        """Create animation widget."""
        # Get the index of the selected mode
        index_str = self.active_modes.value.split(":")[0]
        index = int(index_str) - 1
        # Get the eigenvector of the selected mode
        eigenvector = self.eigenvectors[index]
        # Get the amplitude of the selected mode
        amplitude = self.amplitude.value
        # Get the structure of the selected mode
        structure = self.structure_ase
        # Create a trajectory
        traj = []
        time_array = np.linspace(0.0, 2 * np.pi, 20)
        for time in time_array:
            vibro_atoms = Atoms(
                symbols=structure.symbols,
                positions=structure.positions + amplitude * eigenvector * np.sin(time),
                cell=structure.cell,
                pbc=True,
            )
            supercell = vibro_atoms.repeat(
                (
                    self._supercell[0].value,
                    self._supercell[1].value,
                    self._supercell[2].value,
                )
            )
            traj.append(supercell)
        # Create the animation widget
        self.ngl.clear()
        self.ngl = nv.show_asetraj(traj)
        self.ngl.add_unitcell()
        # Make the animation to be set in a loop
        self.ngl._iplayer.children[0]._playing = False
        self.ngl._iplayer.children[0]._playing = True
        self.ngl._iplayer.children[0]._repeat = True
        self.ngl._set_size("400px", "400px")
        self.ngl.center()
