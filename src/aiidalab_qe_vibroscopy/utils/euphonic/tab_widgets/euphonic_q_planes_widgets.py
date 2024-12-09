import ipywidgets as ipw
import base64
import json
import numpy as np
import plotly.io as pio

from IPython.display import display

from euphonic.cli.utils import (
    _get_debye_waller,
    _get_energy_bins,
)

from aiidalab_qe_vibroscopy.utils.euphonic.base_widgets.euphonic_base_widgets import (
    StructureFactorSettingsBaseWidget,
    COLORSCALE,
    COLORBAR_DICT,
    StructureFactorBasePlotWidget,
)

from monty.json import jsanitize
import plotly.graph_objects as go

from aiidalab_qe_vibroscopy.utils.euphonic.data.structure_factors import (
    blockPrint,
    enablePrint,
)


def produce_Q_section_modes(
    fc,
    h,
    k,
    Q0=np.array([0, 0, 0]),
    n_h=100,
    n_k=100,
    h_extension=1,
    k_extension=1,
    temperature=0,
):
    from euphonic import ureg

    # see get_Q_section
    # h: array vector
    # k: array vector
    # Q0: "point" in Q-space used to build the portion of plane, using also the two vectors h and k.
    # n_h, n_k: number of points along the two directions. or better, the two vectors.

    def get_Q_section(h, k, Q0, n_h, n_k, h_extension, k_extension):
        # every point in the space is Q=Q0+dv1*h+dv2*k, which

        q_list = []
        h_list = []
        k_list = []

        for dv1 in np.linspace(-h_extension, h_extension, n_h):
            for dv2 in np.linspace(-k_extension, k_extension, n_k):
                p = Q0 + dv1 * h + dv2 * k
                q_list.append(p)  # Q list
                h_list.append(dv1)  # *h[0])
                k_list.append(dv2)  # *k[1])

        return np.array(q_list), np.array(h_list), np.array(k_list)

    q_array, h_array, k_array = get_Q_section(
        h, k, Q0, n_h + 1, n_k + 1, h_extension, k_extension
    )

    modes = fc.calculate_qpoint_phonon_modes(qpts=q_array, reduce_qpts=False)

    if temperature > 0:
        blockPrint()
        dw = _get_debye_waller(
            temperature * ureg("K"),
            fc,
            # grid_spacing=(args.grid_spacing * recip_length_unit),
            # **calc_modes_kwargs,
        )
        enablePrint()
    else:
        dw = None

    labels = {
        "q": f"Q0={[np.round(i,3) for i in Q0]}",
        "h": f"h={[np.round(i,3) for i in h]}",
        "k": f"k={[np.round(i,3) for i in k]}",
    }

    return modes, q_array, h_array, k_array, labels, dw


def produce_Q_section_spectrum(
    modes,
    q_array,
    h_array,
    k_array,
    ecenter,
    deltaE=0.5,
    bins=10,
    spectrum_type="coherent",
    dw=None,
    labels=None,
):
    from aiidalab_qe_vibroscopy.utils.euphonic.data.structure_factors import (
        blockPrint,
        enablePrint,
    )

    # bins = 10 # hard coded beacuse here it does not change anything.
    ebins = _get_energy_bins(
        modes, bins + 1, emin=ecenter - deltaE, emax=ecenter + deltaE
    )

    blockPrint()
    if (
        spectrum_type == "coherent"
    ):  # Temperature?? For now let's drop it otherwise it is complicated.
        spectrum = modes.calculate_structure_factor(dw=dw).calculate_sqw_map(ebins)
    elif spectrum_type == "dos":
        spectrum = modes.calculate_dos_map(ebins)

    mu = ecenter
    sigma = (deltaE) / 2

    # Gaussian weights.
    weights = np.exp(-((spectrum.y_data.magnitude - mu) ** 2) / 2 * sigma**2) / np.sqrt(
        2 * np.pi * sigma**2
    )
    av_spec = np.average(spectrum.z_data.magnitude, axis=1, weights=weights[:-1])
    enablePrint()

    return av_spec, q_array, h_array, k_array, labels


class QSectionPlotWidget(StructureFactorBasePlotWidget):
    def render(
        self,
    ):
        if self.rendered:
            return

        super().render()

        self.fig.add_trace(
            go.Heatmap(
                z=self._model.av_spec,
                x=self._model.h_array,
                y=self._model.k_array,
                colorbar=COLORBAR_DICT,
                colorscale=COLORSCALE,  # imported from euphonic_base_widgets
            )
        )

        self.fig.update_layout(
            height=625,
            width=650,
            margin=dict(l=15, r=15, t=15, b=15),
        )

        # self.children.insert(0, ipw.HTML(labels["q"]))
        self.E_units_button.layout.display = "none"
        # self.fig.update_layout(title=labels["q"])
        self.fig["layout"]["xaxis"].update(
            range=[np.min(self._model.h_array), np.max(self._model.h_array)],
            title=r"$\alpha \text{ in } \vec{Q_0} + \alpha \vec{h}$",
        )
        self.fig["layout"]["yaxis"].update(
            range=[np.min(self._model.k_array), np.max(self._model.k_array)],
            title=r"$\beta \text{ in } \vec{Q_0} + \beta \vec{k}$",
        )

    def _update_plot(self):
        self.fig.add_trace(
            go.Heatmap(
                z=self._model.av_spec,
                x=self._model.h_array,
                y=self._model.k_array,
                colorbar=COLORBAR_DICT,
                colorscale=COLORSCALE,  # imported from euphonic_base_widgets
            )
        )

        self.fig.update_layout(
            height=625,
            width=650,
            margin=dict(l=15, r=15, t=15, b=15),
        )
        # self.fig.update_layout(title=labels["q"])
        self.fig["layout"]["xaxis"].update(
            range=[np.min(self._model.h_array), np.max(self._model.h_array)],
            title=r"$\alpha \text{ in } \vec{Q_0} + \alpha \vec{h}$",
        )
        self.fig["layout"]["yaxis"].update(
            range=[np.min(self._model.k_array), np.max(self._model.k_array)],
            title=r"$\beta \text{ in } \vec{Q_0} + \beta \vec{k}$",
        )

        self.fig.data = [self.fig.data[1]]


class QSectionSettingsWidget(StructureFactorSettingsBaseWidget):
    def render(self):
        if self.rendered:
            return

        super().render()

        self.ecenter = ipw.FloatText(
            value=0,
            description="E (meV)",
        )
        ipw.link(
            (self._model, "center_e"),
            (self.ecenter, "value"),
        )
        self.ecenter.observe(self._on_setting_changed, names="value")

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
                child.observe(self._on_setting_changed, names="value")
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

        self.energy_broadening = ipw.FloatText(
            value=0.5,
            description="&Delta;E (meV)",
            tooltip="Energy window in eV",
        )
        ipw.link(
            (self._model, "energy_broadening"),
            (self.energy_broadening, "value"),
        )
        self.energy_broadening.observe(self._on_setting_changed, names="value")

        self.plot_button.disabled = False
        self.plot_button.description = "Plot"
        # self.reset_button.disabled = True
        self.download_button.disabled = True

        # Please note: if you change the order of the widgets below, it will
        # affect the usage of the children[0] below in the full widget.
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
                            # self.specification_intensity,
                            self.ecenter,
                            self.energy_broadening,
                            # self.int_energy_bins, # does not change anything here.
                            self.temperature,
                            self.weight_button,
                        ],
                        layout=ipw.Layout(
                            width="50%",
                        ),
                    ),
                    ipw.VBox(
                        [
                            self.plane_description_widget,
                            self.Q0_widget,
                            self.h_widget,
                            self.k_widget,
                        ],
                        layout=ipw.Layout(
                            width="50%",
                        ),
                    ),
                ],
            ),
        ]

    def _on_vector_changed(self, change=None):
        """
        Update the model.
        """
        self._model.Q0_vec = [i.value for i in self.Q0_vec.children[:-2]]
        self._model.h_vec = [i.value for i in self.h_vec.children]
        self._model.k_vec = [i.value for i in self.k_vec.children]


class QSectionFullWidget(ipw.VBox):
    """
    The Widget to display specifically the Intensity map of Dynamical structure
    factor for single crystal.

    The scattering lengths used in the `produce_bands_weigthed_data` function
    are tabulated (Euphonic/euphonic/data/sears-1992.json)
    and are from Sears (1992) Neutron News 3(3) pp26--37.
    """

    def __init__(self, model):
        self._model = model
        self.rendered = False
        super().__init__()

    def render(self):
        if self.rendered:
            return

        self.title_intensity = ipw.HTML(
            "<h3>Neutron dynamic structure factor in a Q-section visualization</h3>"
        )

        self.settings_intensity = QSectionSettingsWidget(model=self._model)

        # rendering the widget.
        self.settings_intensity.render()

        # specific for the q section
        self.settings_intensity.plot_button.on_click(self._on_plot_button_clicked)
        self.settings_intensity.download_button.on_click(self.download_data)

        self.children = [
            self.title_intensity,
            # self.map_widget,
            self.settings_intensity,
        ]

        self.rendered = True

    def _on_plot_button_clicked(self, change=None):
        self._model._update_qsection_spectra()

        if not hasattr(self, "map_widget"):
            self.map_widget = QSectionPlotWidget(self._model)  # CHANGED
            self.map_widget.render()
            self.children = [
                self.title_intensity,
                self.map_widget,
                self.settings_intensity,
            ]

        self.map_widget._update_plot()

    def download_data(self, _=None):
        """
        Download both the ForceConstants and the spectra json files.
        """
        # force_constants_dict = self.fc.to_dict()

        filename = "Q_section.json"
        my_dict = {}

        my_dict["data"] = self.children[1].fig.data[0].to_plotly_json()

        my_dict.update(
            {
                "h": np.array(
                    [i.value for i in self.settings_intensity.h_vec.children[:-2]]
                ),
                "k": np.array(
                    [i.value for i in self.settings_intensity.k_vec.children[:-2]]
                ),
                "Q0": np.array(
                    [i.value for i in self.settings_intensity.Q0_vec.children[:-2]]
                ),
                "weight": self.settings_intensity.weight_button.value,
                "ecenter": self.settings_intensity.ecenter.value,
                "deltaE": self.settings_intensity.energy_broadening.value,
                "temperature": self.settings_intensity.temperature.value,
            }
        )
        for k in ["h", "k", "Q0", "weight", "ecenter", "deltaE", "temperature"]:
            filename += "_" + k + "_" + str(my_dict[k])

        """
        # FC download:
        json_str = json.dumps(jsanitize(force_constants_dict))
        b64_str = base64.b64encode(json_str.encode()).decode()
        self._download(payload=b64_str, filename="force_constants.json")
        """

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
