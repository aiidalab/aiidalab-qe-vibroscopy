class EuphonicStructureFactorWidget(ipw.VBox):
    """Description.
    
    Collects all the button and widget used to define settings for Neutron dynamic structure factor,
    in all the three cases: single crystal, powder, and q-section....
    """

    def __init__(self, model, spectrum_type = "single_crystal", detached_app = False, **kwargs):
        super().__init__()
        node
        self._model = model
        self._model.spectrum_type = spectrum_type
        self._model.detached_app = detached_app
        self.rendered = False
        
    def render(self):
        """Render the widget.

        This means render the plot button.
        """
        if self.rendered:
            return
        
        self.tab_widget = ipw.Tab()
        self.tab_widget.layout.display = "none"
        self.tab_widget.set_title(0, "Single crystal")
        self.tab_widget.set_title(1, "Powder sample")
        self.tab_widget.set_title(2, "Q-plane view")
        self.tab_widget.children = ()

        self.plot_button = ipw.Button(
            description="Initialise INS data",
            icon="pencil",
            button_style="primary",
            disabled=True,
            layout=ipw.Layout(width="auto"),
        )
        self.plot_button.on_click(self._render_for_real)

        self.loading_widget = LoadingWidget("Loading INS data")
        self.loading_widget.layout.display = "none"

        if not self._model.detached_app:
            self.plot_button.disabled = False
        else:
            self.upload_widget = UploadPhonopyWidget()
            self.upload_widget.reset_uploads.on_click(self._on_reset_uploads_button_clicked)
            self.upload_widget.children[0].observe(self._on_upload_yaml, "_counter")
            self.upload_widget.children[1].observe(self._on_upload_hdf5, "_counter")
            self.children += (self.upload_widget,)
            
        self.download_widget = DownloadYamlHdf5Widget(model=self._model)
        self.download_widget.layout.display = "none"

        self.children += (
            self.plot_button,
            self.loading_widget,
            self.tab_widget,
            self.download_widget,
            )
        
        # NOTE: we initialise here the figure widget, but we do not plot anything yet.
        # this is useful to call the init_view method, which contains the update for the figure.
        self.fig = go.FigureWidget()
        
        self.rendered = True
        
    def _render_for_real(self, change=None):

        self.plot_button.layout.display = "none"
        self.loading_widget.layout.display = "block"
        
        self._init_view()
        
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
        slider_intensity.observe(self._update_intensity_filter, "value")
        specification_intensity = ipw.HTML(
            "(Intensity is relative to the maximum intensity at T=0K)"
        )

        E_units_button = ipw.ToggleButtons(
            options=[
                ("meV", "meV"),
                ("THz", "THz"),
                # ("cm<sup>-1</sup>", "cm-1"),
            ],
            value="meV",
            description="Energy units:",
            disabled=False,
            layout=ipw.Layout(
                width="auto",
            ),
        )
        E_units_button.observe(self._update_energy_units, "value")
        # MAYBE WE LINK ALSO THIS TO THE MODEL? so we can download the data with the preferred units.

        q_spacing = ipw.FloatText(
            value=self._model.q_spacing,
            step=0.001,
            description="q step (1/A)",
            tooltip="q spacing in 1/A",
        )
        ipw.link(
            (self._model, "q_spacing"),
            (self.q_spacing, "value"),
        )
        q_spacing.observe(self._on_setting_change, names="value")

        energy_broadening = ipw.FloatText(
            value=self._model.energy_broadening,
            step=0.01,
            description="&Delta;E (meV)",
            tooltip="Energy broadening in meV",
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
        )
        ipw.link(
            (self._model, "energy_bins"),
            (energy_bins, "value"),
        )
        energy_bins.observe(self._on_setting_change, names="value")

        temperature = ipw.FloatText(
            value=self._model.temperature,
            step=0.01,
            description="T (K)",
            disabled=False,
        )
        ipw.link(
            (self._model, "temperature"),
            (temperature, "value"),
        )
        temperature.observe(self._on_setting_change, names="value")

        weight_button = ipw.ToggleButtons(
            options=[
                ("Coherent", "coherent"),
                ("DOS", "dos"),
            ],
            value=self._model.weighting,
            description="weight:",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "weighting"),
            (weight_button, "value"),
        )
        weight_button.observe(self._on_weight_button_change, names="value")

        plot_button = ipw.Button(
            description="Replot",
            icon="pencil",
            button_style="primary",
            disabled=True,
            layout=ipw.Layout(width="auto"),
        )
        plot_button.observe(self._on_plot_button_change, names="disabled")

        reset_button = ipw.Button(
            description="Reset",
            icon="recycle",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(width="auto"),
        )
        reset_button.on_click(self._reset_settings)

        download_button = ipw.Button(
            description="Download Data and Plot",
            icon="download",
            button_style="primary",
            disabled=False,  # Large files...
            layout=ipw.Layout(width="auto"),
        )
        download_button.on_click(self._download_data)
        
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
            self.custom_kpath_text.observe(self._on_setting_changed, names="value")
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
            self.qmin.observe(self._on_setting_changed, names="value")

            self.qmax = ipw.FloatText(
                step=0.01,
                value=1,
                description="|q|<sub>max</sub> (1/A)",
            )
            ipw.link(
                (self._model, "q_max"),
                (self.qmax, "value"),
            )
            self.qmax.observe(self._on_setting_changed, names="value")

            self.int_npts = ipw.IntText(
                value=100,
                description="npts",
                tooltip="Number of points to be used in the average sphere.",
            )
            ipw.link(
                (self._model, "npts"),
                (self.int_npts, "value"),
            )
            self.int_npts.observe(self._on_setting_changed, names="value")
        # fi self._model.spectrum_type == "powder"
        elif self._model.spectrum_type == "q_planes":
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
        # fi self._model.spectrum_type == "q_planes"
        
        self.children += (
            ...
        )
    
    def _init_view(self, _=None):
        self._model.fetch_data()
        self._update_plot()
        
    def _on_plot_button_change(self, change):
        self.download_button.disabled = not change["new"]

    def _on_weight_button_change(self, change):
        self._model.temperature = 0
        self.temperature.disabled = True if change["new"] == "dos" else False
        self.plot_button.disabled = False

    def _on_setting_change(
        self, change
    ):  # think if we want to do something more evident...
        self.plot_button.disabled = False
        
    def _update_plot(self):
        # update the spectra, i.e. the data contained in the _model.
        # TODO: we need to treat differently the update of intensity and units.
        # they anyway need to modify the data, but no additional spectra re-generation is really needed.
        # so the update_spectra need some more logic, or we call another method.
        self._model.get_spectra()
        
        if not self.rendered:
            # First time we render, we set several layout settings.
            # Layout settings
            self.fig["layout"]["xaxis"].update(
                title=self._model.xlabel,
                range=[min(self._model.x), max(self._model.x)],
            )
            self.fig["layout"]["yaxis"].update(
                title=self._model.ylabel,
                range=[min(self._model.y), max(self._model.y)],
            )
            
            if self.fig.layout.images:
                for image in self.fig.layout.images:
                    image["scl"] = 2  # Set the scale for each image
                    
            self.fig.update_layout(
            height=500,
            width=700,
            margin=dict(l=15, r=15, t=15, b=15),
            )
            # Update x-axis and y-axis to enable autoscaling
            self.fig.update_xaxes(autorange=True)
            self.fig.update_yaxes(autorange=True)

            # Update the layout to enable autoscaling
            self.fig.update_layout(autosize=True)


        heatmap_trace = go.Heatmap(
                z=self._model.z,
                y=(self._model.y),
                x=self._model.x,
                colorbar=COLORBAR_DICT,
                colorscale=COLORSCALE,  # imported from euphonic_base_widgets
            )

        # change the path wants also a change in the labels
        if "ticks_positions" in self._model and "ticks_labels" in self._model:
            self.fig.update_layout(
                xaxis=dict(
                    tickmode="array",
                    tickvals=self._model.ticks_positions,
                    ticktext=self._model.ticks_labels,
                )
            )

        # Add colorbar
        colorbar = heatmap_trace.colorbar
        colorbar.x = 1.05  # Move colorbar to the right
        colorbar.y = 0.5  # Center colorbar vertically

        # Add heatmap trace to figure
        self.fig.add_trace(heatmap_trace)
        self.fig.data = [self.fig.data[1]]

    def _reset_settings(self, _):
        self._model.reset()
        
    def _download_data(self, _=None):
        data, filename = self._model.prepare_data_for_download()
        self._download(data, filename)
    
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

    def _on_vector_changed(self, change=None):
        """
        Update the model. Specific to qplanes case.
        """
        self._model.Q0_vec = [i.value for i in self.Q0_vec.children[:-2]]
        self._model.h_vec = [i.value for i in self.h_vec.children]
        self._model.k_vec = [i.value for i in self.k_vec.children]
