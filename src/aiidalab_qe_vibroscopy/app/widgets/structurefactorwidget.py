import ipywidgets as ipw
import numpy as np
import plotly.graph_objs as go

from IPython.display import display

from aiidalab_qe.common.infobox import InfoBox

from aiidalab_qe_vibroscopy.app.widgets.euphonicmodel import EuphonicResultsModel

COLORSCALE = "Viridis"  # we should allow more options
COLORBAR_DICT = dict(orientation="v", showticklabels=False, x=1, thickness=10, len=0.4)


class EuphonicStructureFactorWidget(ipw.VBox):
    """The true Euphonic widget (to not be confused with the collective EuphonicWidget).

    These is a three-fold widget: it can be a single-crystal, a powder or a q-plane widget.
    This is decided by the spectrum_type attribute of the model.
    """

    def __init__(
        self,
        model: EuphonicResultsModel,
        node=None,
        spectrum_type="single_crystal",
        **kwargs,
    ):
        super().__init__()
        self._model = model
        if node:
            self._model.vibro = node
        if not hasattr(self._model, "spectrum_type"):
            self._model.spectrum_type = spectrum_type

        self.rendered = False

    def render(self):
        """Renders the widget.

        Basically, we have three main parts:
        - a toggle button to display plot info, parameters descriptions and so on
        - parameters
        - plot

        we first create the common widgets (for the three spectrum_type cases), then we add the specific ones
        via if-elif-else.
        """
        if self.rendered:
            return

        # 1 - info part
        self.about_toggle = ipw.ToggleButton(
            layout=ipw.Layout(width="20%"),
            button_style="",
            icon="info",
            value=False,
            description="Plot info",
            tooltip="Info on the parameters and plots",
            disabled=False,
        )
        self.about_toggle.observe(self.display_info_legend, names="value")

        self.info_legend = ipw.HTML("")
        ipw.dlink(
            (self._model, "info_legend_text"),
            (self.info_legend, "value"),
        )
        self.info_legend_infobox = InfoBox(
            children=[self.info_legend],
        )
        self.info_legend_infobox.layout.display = "none"
        self._model.generate_info_legend()

        # 2 - parameters part
        slider_intensity = ipw.FloatRangeSlider(
            value=[1, 100],  # Default selected interval
            min=0,
            max=100,
            step=0.1,
            orientation="horizontal",
            readout=True,
            readout_format=".1f",
            layout=ipw.Layout(
                width="400px",
            ),
        )
        ipw.link(
            (slider_intensity, "value"),
            (self._model, "intensity_filter"),
        )
        slider_intensity.observe(self._update_intensity_filter, "value")
        specification_intensity = ipw.HTML("Intensity window (%):")

        E_units_dropdown = ipw.Dropdown(
            options=[
                ("meV", "meV"),
                ("THz", "THz"),
                ("1/cm", "1/cm"),
            ],
            value="meV",
            description="Energy units:",
            disabled=False,
            layout=ipw.Layout(
                width="auto",
            ),
        )
        ipw.link(
            (E_units_dropdown, "value"),
            (self._model, "energy_units"),
        )
        E_units_dropdown.observe(self._update_energy_units, "value")

        q_spacing = ipw.BoundedFloatText(
            value=self._model.q_spacing,
            step=0.01,
            min=0.01,
            description="q step (1/A)",
            tooltip="q spacing in 1/A",
            layout=ipw.Layout(
                width="auto",
            ),
            continuous_update=True,
        )
        ipw.link(
            (self._model, "q_spacing"),
            (q_spacing, "value"),
        )
        q_spacing.observe(self._on_setting_change, names="value")

        self.energy_broadening = ipw.BoundedFloatText(
            value=self._model.energy_broadening,
            step=0.01,
            min=0,
            description="&Delta;E (meV)",
            tooltip="Energy broadening in meV",
            layout=ipw.Layout(
                width="auto",
            ),
            continuous_update=True,
        )
        ipw.link(
            (self._model, "energy_broadening"),
            (self.energy_broadening, "value"),
        )
        self.energy_broadening.observe(self._on_setting_change, names="value")

        ebins = ipw.BoundedIntText(
            value=self._model.ebins,
            min=1,
            max=5000,
            description="#E bins",
            tooltip="Number of energy bins",
            layout=ipw.Layout(
                width="auto",
            ),
            continuous_update=True,
        )
        ipw.link(
            (self._model, "ebins"),
            (ebins, "value"),
        )
        ebins.observe(self._on_setting_change, names="value")

        self.temperature = ipw.BoundedFloatText(
            value=self._model.temperature,
            step=0.01,
            min=0,
            max=1000,
            description="T (K)",
            disabled=False,
            layout=ipw.Layout(
                width="auto",
            ),
            continuous_update=True,
        )
        ipw.link(
            (self._model, "temperature"),
            (self.temperature, "value"),
        )
        self.temperature.observe(self._on_setting_change, names="value")

        weight_button = ipw.Dropdown(
            options=[
                ("S(Q, ω)", "coherent"),
                ("DOS map", "dos"),
            ],
            value=self._model.weighting,
            description="Plot mode:",
            disabled=False,
            style={"description_width": "initial"},
            layout=ipw.Layout(
                width="auto",
            ),
        )
        ipw.link(
            (self._model, "weighting"),
            (weight_button, "value"),
        )
        weight_button.observe(self._on_weight_button_change, names="value")

        self.plot_button = ipw.Button(
            description="Replot",
            icon="pencil",
            button_style="primary",
            disabled=True,
            layout=ipw.Layout(width="auto"),
        )
        self.plot_button.on_click(self._update_plot)

        reset_button = ipw.Button(
            description="Reset",
            icon="recycle",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(width="auto"),
        )
        reset_button.on_click(self._reset_settings)

        self.download_button = ipw.Button(
            description="Download Data",
            icon="download",
            button_style="primary",
            disabled=False,  # Large files...
            layout=ipw.Layout(width="auto"),
        )
        self.download_button.on_click(self._model._download_data)
        ipw.dlink(
            (self.plot_button, "disabled"),
            (self.download_button, "disabled"),
            lambda x: not x,
        )

        self._init_view()  # this will generate the self.figure_container

        self.children += (
            self.about_toggle,
            self.info_legend_infobox,
            ipw.HBox(
                [
                    specification_intensity,
                    slider_intensity,
                ],
                layout=ipw.Layout(
                    justify_content="flex-start",
                    # margin="10px 0",
                ),
            ),
            ipw.HBox(
                [
                    ipw.VBox(
                        [
                            E_units_dropdown,
                            q_spacing,
                            self.energy_broadening,
                            ebins,
                            self.temperature,
                            weight_button,
                            self.plot_button,
                            reset_button,
                            self.download_button,
                        ],
                        layout=ipw.Layout(
                            justify_content="flex-start", max_width="250px"
                        ),
                    ),
                    self.figure_container,
                ],
            ),
        )

        if self._model.spectrum_type == "single_crystal":
            self.custom_kpath_text = ipw.Text(
                value="",
                description="Custom path (rlu):",
                style={"description_width": "initial"},
                continuous_update=True,
            )
            custom_style = '<style>.custom-font { font-family: "Monospaced"; font-size: 16px; }</style>'
            display(ipw.HTML(custom_style))
            self.custom_kpath_text.add_class("custom-font")
            ipw.link(
                (self._model, "custom_kpath"),
                (self.custom_kpath_text, "value"),
            )
            self.custom_kpath_text.observe(self._on_setting_change, names="value")

            self.children += (self.custom_kpath_text,)

        elif self._model.spectrum_type == "powder":
            self.qmin = ipw.BoundedFloatText(
                value=0,
                min=0,
                max=1,
                description="|q|<sub>min</sub> (1/A)",
                continuous_update=True,
            )
            ipw.link(
                (self._model, "q_min"),
                (self.qmin, "value"),
            )
            self.qmin.observe(self._on_setting_change, names="value")

            self.qmax = ipw.BoundedFloatText(
                step=0.01,
                value=1,
                min=0,
                max=1,
                description="|q|<sub>max</sub> (1/A)",
                continuous_update=True,
            )
            ipw.link(
                (self._model, "q_max"),
                (self.qmax, "value"),
            )
            self.qmax.observe(self._on_setting_change, names="value")

            self.children += (
                ipw.HBox(
                    [
                        self.qmin,
                        self.qmax,
                    ],
                ),
            )

        elif self._model.spectrum_type == "q_planes":
            q_spacing.layout.display = "none"

            self.ecenter = ipw.FloatText(
                value=0,
                description="E cut (meV)",
                continuous_update=True,
            )
            ipw.link(
                (self._model, "center_e"),
                (self.ecenter, "value"),
            )
            ipw.dlink(
                (self._model, "energy_units"),
                (self.ecenter, "description"),
                transform=lambda x: f"E ({x})",
            )
            self.ecenter.observe(self._on_setting_change, names="value")

            self.plane_description_widget = ipw.HTML(
                """
                <div style="padding-top: 0px; padding-bottom: 0px; line-height: 140%;">
                    <b>Q-plane definition</b>: <br>
                </div>
                """
            )

            self.Q0_vec = ipw.HBox(
                [
                    ipw.FloatText(
                        value=0,
                        layout={"width": "60px"},
                        continuous_update=True,
                    )
                    for j in range(3)
                ]
                + [
                    ipw.HTML(
                        "N<sup>h</sup><sub>q</sub>, N<sup>k</sup><sub>q</sub> &darr;",
                        layout={"width": "60px"},
                    ),
                    ipw.HTML(r"&alpha;, &beta; &darr;", layout={"width": "60px"}),
                ]
            )

            self.h_vec = ipw.HBox(
                [
                    ipw.FloatText(
                        value=1,
                        layout={"width": "60px"},
                        continuous_update=True,
                    )  # coordinates
                    for j in range(3)
                ]
                + [
                    ipw.IntText(value=100, layout={"width": "60px"}),
                    ipw.IntText(value=1, layout={"width": "60px"}),
                ]  # number of points along this dir, i.e. n_h; and multiplicative factor alpha
            )
            self.k_vec = ipw.HBox(
                [
                    ipw.FloatText(
                        value=1,
                        layout={"width": "60px"},
                        continuous_update=True,
                    )
                    for j in range(3)
                ]
                + [
                    ipw.IntText(value=100, layout={"width": "60px"}),
                    ipw.IntText(value=1, layout={"width": "60px"}),
                ]
            )

            for vec in [self.Q0_vec, self.h_vec, self.k_vec]:
                for child in vec.children:
                    child.observe(self._on_setting_change, names="value")
                    child.observe(self._on_vector_changed, names="value")

            self.Q0_widget = ipw.HBox(
                [ipw.HTML("Q<sub>0</sub>: ", layout={"width": "20px"}), self.Q0_vec]
            )
            self.h_widget = ipw.HBox(
                [ipw.HTML("h:  ", layout={"width": "20px"}), self.h_vec]
            )
            self.k_widget = ipw.HBox(
                [ipw.HTML("k:  ", layout={"width": "20px"}), self.k_vec]
            )

            self.children += (
                self.ecenter,
                self.plane_description_widget,
                self.Q0_widget,
                self.h_widget,
                self.k_widget,
            )

        # RENDERING IS DONE, SO:
        self.rendered = True

    def _init_view(self, _=None):
        # for safety, we fetch the data again (should have happened already in the EuophonicWidget).
        # if already there, this model method will not do anything.
        self._model.fetch_data()

        # we need to initialize the figure container
        if not hasattr(self, "fig"):
            self.fig = go.FigureWidget()
            self.fig.update_layout(
                margin=dict(l=20, r=0, t=0, b=20),
                yaxis=dict(autorange=True),
            )

            self.figure_container = ipw.VBox([self.fig])

        # we populate the plot calling the update method
        self._update_plot()

    def _on_weight_button_change(self, change):
        self._model.temperature = 0
        self.temperature.disabled = True if change["new"] == "dos" else False
        self.plot_button.disabled = False

    def _on_setting_change(
        self, change
    ):  # think if we want to do something more evident...
        self.plot_button.disabled = False

    def _update_plot(self, _=None):
        # update the spectra, i.e. the data to be plotted contained in the _model.
        self._model.get_spectra()

        if self._model.spectrum_type == "q_planes":
            # hide figure until we have the data
            self.figure_container.layout.display = (
                "none" if not self.rendered else "block"
            )
            self.plot_button.disabled = self.rendered

        self.fig.update_layout(yaxis_title=self._model.ylabel)

        # a specific q-path wants the appropriate labels
        if hasattr(self._model, "ticks_positions") and hasattr(
            self._model, "ticks_labels"
        ):
            self.fig.update_layout(
                xaxis=dict(
                    tickmode="array",
                    tickvals=self._model.ticks_positions,
                    ticktext=self._model.ticks_labels,
                )
            )
        elif hasattr(self._model, "xlabel"):
            self.fig.update_layout(xaxis_title=self._model.xlabel)

        # Generate the plot object, in this case a heatmap
        heatmap_trace = go.Heatmap(
            z=self._model.z,
            y=self._model.y,
            x=self._model.x,
            colorscale=COLORSCALE,
        )

        # Add heatmap trace to figure
        self.fig.add_trace(heatmap_trace)
        self.fig.data = [self.fig.data[-1]]

        if self.rendered:
            self._update_intensity_filter()

        self.plot_button.disabled = True

    def _update_intensity_filter(self, change=None):
        # the value of the intensity slider is in fractions of the max.
        # NOTE: we do this here, as we do not want to replot. Reason is that
        # the data will not change! so we don't need to invoke the model.
        self.fig.data[0].zmax = (
            self._model.intensity_filter[1] * np.max(self.fig.data[0].z) / 100
        )  # above this, it is all yellow, i.e. max intensity.
        self.fig.data[0].zmin = (
            self._model.intensity_filter[0] * np.max(self.fig.data[0].z) / 100
        )  # below this, it is all blue, i.e. zero intensity

    def _update_energy_units(self, change):
        """Updating the energy units

        we need to update both the plot and the energy broadening widget.
        Also, we need to avoid to enable the plot button if it was already disabled:
        we are not changing anything apart the units.
        """
        replot_was_off = self.plot_button.disabled

        with self.hold_trait_notifications():
            self._model._update_energy_units(
                old_units=change["old"], new_units=change["new"]
            )

            self.fig.data[0].y = self._model.y
            self.fig.update_layout(yaxis_title=self._model.ylabel)
            self.energy_broadening.description = (
                f"&Delta;E ({self._model.energy_units})"
            )

            # putting off again the replot if it was off.
            self.plot_button.disabled = replot_was_off

    def _reset_settings(self, _):
        self._model.reset()

    def _on_vector_changed(self, change=None):
        """
        Update the model. Specific to qplanes case.
        """
        self._model.Q0_vec = [i.value for i in self.Q0_vec.children[:-2]]
        self._model.h_vec = [i.value for i in self.h_vec.children]
        self._model.k_vec = [i.value for i in self.k_vec.children]

    def display_info_legend(self, change):
        self.info_legend_infobox.layout.display = "block" if change["new"] else "none"
