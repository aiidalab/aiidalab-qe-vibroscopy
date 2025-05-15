import ipywidgets as ipw
from aiidalab_qe_vibroscopy.app.widgets.ramanmodel import RamanModel
import plotly.graph_objects as go
from aiidalab_widgets_base.utils import StatusHTML
from IPython.display import HTML, clear_output, display
from weas_widget import WeasWidget


class RamanWidget(ipw.VBox):
    """
    Widget to display Raman properties Tab
    """

    def __init__(
        self, model: RamanModel, node: None, input_structure, spectrum_type, **kwargs
    ):
        super().__init__(
            children=[ipw.HTML("Loading Raman data...")],
            **kwargs,
        )
        self._model = model
        self._model.spectrum_type = spectrum_type
        self._model.vibro = node
        self._model.input_structure = input_structure
        self._model.fetch_data()
        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.guiConfig = {
            "components": {"enabled": True, "atomsControl": True, "buttons": True},
            "buttons": {
                "enabled": True,
                "fullscreen": True,
                "download": True,
                "measurement": True,
            },
        }

        self.plot_type = ipw.ToggleButtons(
            description="Spectrum type:",
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self._model, "plot_type_options"),
            (self.plot_type, "options"),
        )
        ipw.link(
            (self._model, "plot_type"),
            (self.plot_type, "value"),
        )
        self.plot_type.observe(self._on_plot_type_change, names="value")
        self.temperature = ipw.FloatText(
            description="Temperature (K):",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "temperature"),
            (self.temperature, "value"),
        )
        self.frequency_laser = ipw.FloatText(
            description="Laser frequency (nm):",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "frequency_laser"),
            (self.frequency_laser, "value"),
        )
        self.pol_incoming = ipw.Text(
            description="Incoming polarization:",
            style={"description_width": "initial"},
            layout=ipw.Layout(visibility="hidden"),
        )
        ipw.link(
            (self._model, "pol_incoming"),
            (self.pol_incoming, "value"),
        )
        self.pol_outgoing = ipw.Text(
            description="Outgoing polarization:",
            style={"description_width": "initial"},
            layout=ipw.Layout(visibility="hidden"),
        )
        ipw.link(
            (self._model, "pol_outgoing"),
            (self.pol_outgoing, "value"),
        )
        self.plot_button = ipw.Button(
            description="Update Plot",
            icon="pencil",
            button_style="primary",
            layout=ipw.Layout(width="auto"),
        )
        self.use_nac_direction = ipw.Checkbox(
            description="Use non-analytical direction(NAC)",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "use_nac_direction"),
            (self.use_nac_direction, "value"),
        )
        self.use_nac_direction.observe(self._on_nac_direction_change, names="value")
        self.help_nac_direction = ipw.HTML(
            value="""<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
                The NAC direction should match the light propagation direction, which is perpendicular to the polarization direction, and it should be defined in Cartesian coordinates.
                </div>"""
        )
        self.nac_direction = ipw.Text(
            description="NAC direction:",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "nac_direction"),
            (self.nac_direction, "value"),
        )
        self.plot_button.on_click(self._on_plot_button_click)
        self.download_button = ipw.Button(
            description="Download Data",
            icon="download",
            button_style="primary",
            layout=ipw.Layout(width="auto"),
        )
        self.download_button.on_click(self._model.download_data)
        self._wrong_syntax = StatusHTML(clear_after=8)

        self.broadening = ipw.FloatText(
            description="Broadening (cm<sup>-1</sup>):",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "broadening"),
            (self.broadening, "value"),
        )

        self.separate_polarizations = ipw.Checkbox(
            description="Separate polarized and depolarized intensities",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "separate_polarizations"),
            (self.separate_polarizations, "value"),
        )
        self.spectrum = go.FigureWidget(
            layout=go.Layout(
                title=dict(text="Powder Raman spectrum"),
                barmode="overlay",
                xaxis=dict(
                    title="Wavenumber (cm<sup>-1</sup>)",
                    nticks=0,
                ),
                yaxis=dict(
                    title="Intensity (arb. units)",
                ),
                height=500,
                width=700,
                plot_bgcolor="white",
            )
        )

        # Active Modes
        self.modes_table = ipw.Output()
        self.animation = ipw.Output()

        self.active_modes = ipw.Dropdown(
            description="Select mode:",
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self._model, "active_modes_options"),
            (self.active_modes, "options"),
        )
        ipw.link(
            (self._model, "active_mode"),
            (self.active_modes, "value"),
        )
        self.amplitude = ipw.FloatText(
            description="Amplitude :",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "amplitude"),
            (self.amplitude, "value"),
        )
        self._supercell = [
            ipw.BoundedIntText(min=1, layout={"width": "40px"}),
            ipw.BoundedIntText(min=1, layout={"width": "40px"}),
            ipw.BoundedIntText(min=1, layout={"width": "40px"}),
        ]
        for i, widget in enumerate(self._supercell):
            ipw.link(
                (self._model, f"supercell_{i}"),
                (widget, "value"),
            )

        self.supercell_selector = ipw.HBox(
            [
                ipw.HTML(
                    description="Super cell:", style={"description_width": "initial"}
                )
            ]
            + self._supercell
        )
        # WeasWidget Setting
        self.weas = WeasWidget(
            guiConfig=self.guiConfig, viewerStyle={"width": "800px", "height": "600px"}
        )
        self.weas.from_ase(self._model.input_structure)
        self.weas.avr.model_style = 1
        self.weas.avr.color_type = "JMOL"

        widget_list = [
            self.active_modes,
            self.amplitude,
            self._supercell[0],
            self._supercell[1],
            self._supercell[2],
        ]
        for elem in widget_list:
            elem.observe(self._select_active_mode, names="value")

        self.children = [
            ipw.HTML(f"<h3>{self._model.spectrum_type} spectroscopy</h3>"),
            ipw.HTML(
                """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
                Select the type spectrum to plot.
                </div>"""
            ),
            self.plot_type,
            self.temperature,
            self.frequency_laser,
            self.broadening,
            self.separate_polarizations,
            self.use_nac_direction,
            self.help_nac_direction,
            self.nac_direction,
            self.pol_incoming,
            self.pol_outgoing,
            self._wrong_syntax,
            ipw.HBox([self.plot_button, self.download_button]),
            self.spectrum,
            ipw.HBox(
                [
                    ipw.VBox(
                        [
                            ipw.HTML(
                                value=f"<b>{self._model.spectrum_type} active modes</b>"
                            ),
                            self.modes_table,
                        ]
                    ),
                    ipw.VBox(
                        [
                            self.active_modes,
                            self.amplitude,
                            self.supercell_selector,
                        ],
                    ),
                ]
            ),
            self.animation,
        ]

        self.rendered = True

        self._initial_view()

    def _initial_view(self):
        self._model.update_data()
        if self._model.spectrum_type == "IR":
            self.temperature.layout.display = "none"
            self.frequency_laser.layout.display = "none"
            self.pol_outgoing.layout.display == "none"
            self.separate_polarizations.layout.display = "none"

        self.nac_direction.layout.display = "none"
        self.help_nac_direction.layout.display = "none"

        self.spectrum.add_scatter(
            x=self._model.frequencies, y=self._model.intensities, name=""
        )
        self.spectrum.layout.title.text = f"Powder {self._model.spectrum_type} spectrum"
        self.modes_table.layout = {
            "overflow": "auto",
            "height": "200px",
            "width": "150px",
        }
        self.weas = self._model.set_vibrational_mode_animation(self.weas)

        with self.animation:
            clear_output()
            display(self.weas)

        with self.modes_table:
            clear_output()
            display(HTML(self._model.modes_table()))

    def _on_nac_direction_change(self, change):
        if change["new"]:
            self.nac_direction.layout.display = "block"
            self.help_nac_direction.layout.display = "block"
        else:
            self.nac_direction.layout.display = "none"
            self.help_nac_direction.layout.display = "none"

    def _on_plot_type_change(self, change):
        if change["new"] == "single_crystal":
            self.pol_incoming.layout.visibility = "visible"
            if self._model.spectrum_type == "Raman":
                self.pol_outgoing.layout.visibility = "visible"
            else:
                self.separate_polarizations.layout.display = "none"
            self.separate_polarizations.layout.visibility = "hidden"
        else:
            self.pol_incoming.layout.visibility = "hidden"
            self.pol_outgoing.layout.visibility = "hidden"
            self.separate_polarizations.layout.visibility = "visible"

    def _on_plot_button_click(self, _):
        _, incoming_syntax_ok = self._model._check_inputs_correct(
            self.pol_incoming.value
        )
        _, outgoing_syntax_ok = self._model._check_inputs_correct(
            self.pol_outgoing.value
        )
        _, nac_syntax_ok = self._model._check_inputs_correct(self.nac_direction.value)
        if not (incoming_syntax_ok and outgoing_syntax_ok and nac_syntax_ok):
            self._wrong_syntax.message = """
                <div class='alert alert-danger'>
                    ERROR: Invalid syntax for polarization directions.
                </div>
            """
            return
        self._model.update_data()
        self._model.update_plot(self.spectrum)

    def _select_active_mode(self, _):
        self.weas = self._model.set_vibrational_mode_animation(self.weas)
        with self.animation:
            clear_output()
            display(self.weas)
