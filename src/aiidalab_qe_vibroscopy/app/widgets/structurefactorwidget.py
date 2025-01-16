import ipywidgets as ipw
from IPython.display import display

import numpy as np

# ADD ALL THE IMPORTS.
import plotly.graph_objs as go

from aiidalab_qe_vibroscopy.app.widgets.euphonicmodel import EuphonicResultsModel

COLORSCALE = "Viridis"  # we should allow more options
COLORBAR_DICT = dict(orientation="v", showticklabels=False, x=1, thickness=10, len=0.4)


class EuphonicStructureFactorWidget(ipw.VBox):
    """Description.

    Collects all the button and widget used to define settings for Neutron dynamic structure factor,
    in all the three cases: single crystal, powder, and q-section....
    """

    def __init__(
        self,
        model: EuphonicResultsModel,
        node=None,
        spectrum_type="single_crystal",
        detached_app=False,
        **kwargs,
    ):
        super().__init__()
        self._model = model
        if node:
            self._model.vibro = node
        self._model.spectrum_type = spectrum_type
        self._model.detached_app = detached_app
        self.rendered = False

    def render(self):
        """Render the widget.

        This means render the plot button.
        """
        if self.rendered:
            return

        slider_intensity = ipw.FloatRangeSlider(
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
        ipw.link(
            (slider_intensity, "value"),
            (self._model, "intensity_filter"),
        )
        slider_intensity.observe(self._update_intensity_filter, "value")

        specification_intensity = ipw.HTML(
            "Intensity window (relative to the maximum intensity at T=0K):"
        )

        E_units_ddown = ipw.Dropdown(
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
            (E_units_ddown, "value"),
            (self._model, "energy_units"),
        )
        E_units_ddown.observe(self._update_energy_units, "value")

        q_spacing = ipw.FloatText(
            value=self._model.q_spacing,
            step=0.001,
            description="q step (1/A)",
            tooltip="q spacing in 1/A",
            layout=ipw.Layout(
                width="auto",
            ),
        )
        ipw.link(
            (self._model, "q_spacing"),
            (q_spacing, "value"),
        )
        q_spacing.observe(self._on_setting_change, names="value")

        energy_broadening = ipw.FloatText(
            value=self._model.energy_broadening,
            step=0.01,
            description="&Delta;E (meV)",
            tooltip="Energy broadening in meV",
            layout=ipw.Layout(
                width="auto",
            ),
        )
        ipw.link(
            (self._model, "energy_broadening"),
            (energy_broadening, "value"),
        )
        energy_broadening.observe(self._on_setting_change, names="value")

        energy_bins = ipw.IntText(
            value=self._model.energy_bins,
            description="#E bins",
            tooltip="Number of energy bins",
            layout=ipw.Layout(
                width="auto",
            ),
        )
        ipw.link(
            (self._model, "energy_bins"),
            (energy_bins, "value"),
        )
        energy_bins.observe(self._on_setting_change, names="value")

        self.temperature = ipw.FloatText(
            value=self._model.temperature,
            step=0.01,
            description="T (K)",
            disabled=False,
            layout=ipw.Layout(
                width="auto",
            ),
        )
        ipw.link(
            (self._model, "temperature"),
            (self.temperature, "value"),
        )
        self.temperature.observe(self._on_setting_change, names="value")

        weight_button = ipw.ToggleButtons(
            options=[
                ("Coherent", "coherent"),
                ("DOS", "dos"),
            ],
            value=self._model.weighting,
            description="weight:",
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

        self._init_view()  # generate the self.figure_container

        self.children += (
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
                            E_units_ddown,
                            q_spacing,
                            energy_broadening,
                            energy_bins,
                            self.temperature,
                            weight_button,
                            self.plot_button,
                            reset_button,
                            self.download_button,
                        ],
                        layout=ipw.Layout(
                            justify_content="flex-start", max_width="200px"
                        ),
                    ),
                    self.figure_container,
                ],
            ),
        )

        if self._model.spectrum_type == "single_crystal":
            self.custom_kpath_description = ipw.HTML(
                """
            <div style="padding-top: 0px; padding-bottom: 0px; line-height: 140%;">
                <b>Custom q-points path for the structure factor</b>: <br>
                we can provide it via a specific format: <br>
                (1) each linear path should be divided by '|'; <br>
                (2) each path is composed of 'qxi qyi qzi - qxf qyf qzf' where qxi and qxf are, respectively,
                the start and end x-coordinate of the q direction, in reciprocal lattice units (rlu).<br>
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
            ipw.link(
                (self._model, "custom_kpath"),
                (self.custom_kpath_text, "value"),
            )
            self.custom_kpath_text.observe(self._on_setting_change, names="value")

            self.children += (
                self.custom_kpath_description,
                self.custom_kpath_text,
            )
        # fi self._model.spectrum_type == "single_crystal"
        elif self._model.spectrum_type == "powder":
            self.qmin = ipw.FloatText(
                value=0,
                description="|q|<sub>min</sub> (1/A)",
            )
            ipw.link(
                (self._model, "q_min"),
                (self.qmin, "value"),
            )
            self.qmin.observe(self._on_setting_change, names="value")

            self.qmax = ipw.FloatText(
                step=0.01,
                value=1,
                description="|q|<sub>max</sub> (1/A)",
            )
            ipw.link(
                (self._model, "q_max"),
                (self.qmax, "value"),
            )
            self.qmax.observe(self._on_setting_change, names="value")

            # int_npts fixed to 500 for now.
            # self.int_npts = ipw.IntText(
            #    value=100,
            #    description="npts",
            #    tooltip="Number of points to be used in the average sphere.",
            # )
            # ipw.link(
            #    (self._model, "npts"),
            #    (self.int_npts, "value"),
            # )
            # self.int_npts.observe(self._on_setting_change, names="value")
            self.children += (
                ipw.HBox(
                    [
                        self.qmin,
                        self.qmax,
                    ],
                ),
                # self.int_npts,
            )

        # fi self._model.spectrum_type == "powder"
        elif self._model.spectrum_type == "q_planes":
            E_units_ddown.layout.display = "none"
            q_spacing.layout.display = "none"

            self.ecenter = ipw.FloatText(
                value=0,
                description="E (meV)",
            )
            ipw.link(
                (self._model, "center_e"),
                (self.ecenter, "value"),
            )
            self.ecenter.observe(self._on_setting_change, names="value")

            self.plane_description_widget = ipw.HTML(
                """
                <div style="padding-top: 0px; padding-bottom: 0px; line-height: 140%;">
                    <b>Q-plane definition</b>: <br>
                    To define a plane in the reciprocal space, <br>
                    you should define a point in the reciprocal space, Q<sub>0</sub>,
                    and two vectors h&#8407; and k&#8407;. Then, each Q point is defined as: Q = Q<sub>0</sub> + &alpha;*h&#8407; + &beta;*k&#8407;. <br>
                    Then you can select the number of q points in both directions and the &alpha; and &beta; parameters. <br>
                    Coordinates are reciprocal lattice units (rlu).
                </div>
                """
            )

            self.Q0_vec = ipw.HBox(
                [ipw.FloatText(value=0, layout={"width": "60px"}) for j in range(3)]
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
                    ipw.FloatText(value=1, layout={"width": "60px"})  # coordinates
                    for j in range(3)
                ]
                + [
                    ipw.IntText(value=100, layout={"width": "60px"}),
                    ipw.IntText(value=1, layout={"width": "60px"}),
                ]  # number of points along this dir, i.e. n_h; and multiplicative factor alpha
            )
            self.k_vec = ipw.HBox(
                [ipw.FloatText(value=1, layout={"width": "60px"}) for j in range(3)]
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
        # fi self._model.spectrum_type == "q_planes"

        self.rendered = True

    def _init_view(self, _=None):
        self._model.fetch_data()
        if not hasattr(self, "fig"):
            self.fig = go.FigureWidget()
            self.fig.update_layout(margin=dict(l=20, r=0, t=0, b=20))
            self.figure_container = ipw.VBox([self.fig])
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
        # update the spectra, i.e. the data contained in the _model.

        self._model.get_spectra()

        if self._model.spectrum_type == "q_planes":
            # hide figure until we have the data
            self.figure_container.layout.display = (
                "none" if not self.rendered else "block"
            )
            self.plot_button.disabled = self.rendered

        self.fig.update_layout(yaxis_title=self._model.ylabel)

        # changing the path wants also a change in the labels
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

        heatmap_trace = go.Heatmap(
            z=self._model.z,
            y=self._model.y,
            x=self._model.x,
            # colorbar=COLORBAR_DICT,
            colorscale=COLORSCALE,
        )
        # Add colorbar
        # Do we need it?
        # colorbar = heatmap_trace.colorbar
        # colorbar.x = 1.05  # Move colorbar to the right
        # colorbar.y = 0.5  # Center colorbar vertically

        # Add heatmap trace to figure
        self.fig.add_trace(heatmap_trace)
        self.fig.data = [self.fig.data[-1]]

        if self.rendered:
            self._update_intensity_filter()

        self.plot_button.disabled = True

    def _update_intensity_filter(self, change=None):
        # the value of the intensity slider is in fractions of the max.
        # NOTE: we do this here, as we do not want to replot.
        self.fig.data[0].zmax = (
            self._model.intensity_filter[1] * np.max(self.fig.data[0].z) / 100
        )  # above this, it is all yellow, i.e. max intensity.
        self.fig.data[0].zmin = (
            self._model.intensity_filter[0] * np.max(self.fig.data[0].z) / 100
        )  # below this, it is all blue, i.e. zero intensity

    def _update_energy_units(self, change):
        # the value of the intensity slider is in fractions of the max.
        self._model._update_energy_units(
            new=change["new"],
            old=change["old"],
        )
        # Update x-axis and y-axis to enable autoscaling

        self.fig.data[0].y = self._model.y
        self.fig.update_layout(yaxis_title=self._model.ylabel)

    def _reset_settings(self, _):
        self._model.reset()

    def _on_vector_changed(self, change=None):
        """
        Update the model. Specific to qplanes case.
        """
        self._model.Q0_vec = [i.value for i in self.Q0_vec.children[:-2]]
        self._model.h_vec = [i.value for i in self.h_vec.children]
        self._model.k_vec = [i.value for i in self.k_vec.children]
