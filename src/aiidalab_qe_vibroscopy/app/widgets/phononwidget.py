import ipywidgets as ipw

from aiidalab_qe.common.widgets import LoadingWidget
from aiidalab_qe_vibroscopy.app.widgets.phononmodel import PhononModel

import plotly.graph_objects as go
from aiidalab_qe.common.bands_pdos.bandpdosplotly import BandsPdosPlotly


class PhononWidget(ipw.VBox):
    """
    Widget for displaying phonon properties results
    """

    def __init__(self, model: PhononModel, node: None, **kwargs):
        super().__init__(
            children=[LoadingWidget("Loading widgets")],
            **kwargs,
        )
        self._model = model
        self._model.vibro = node
        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.bandspdos_download_button = ipw.Button(
            description="Download phonon bands and dos data",
            icon="pencil",
            button_style="primary",
            layout=ipw.Layout(width="300px"),
        )
        self.bandspdos_download_button.on_click(self._model.download_bandspdos_data)

        self.thermal_plot = go.FigureWidget(
            layout=go.Layout(
                title=dict(text="Thermal properties"),
                barmode="overlay",
            )
        )

        self.thermo_download_button = ipw.Button(
            description="Download thermal properties data",
            icon="pencil",
            button_style="primary",
            layout=ipw.Layout(width="300px"),
        )
        self.thermo_download_button.on_click(self._model.download_thermo_data)

        self.children = [
            self.bandspdos_download_button,
            self.thermal_plot,
            self.thermo_download_button,
        ]

        self.rendered = True
        self._init_view()

    def _init_view(self):
        self._model.fetch_data()
        self.bands_pdos = BandsPdosPlotly(
            bands_data=self._model.bands_data, pdos_data=self._model.pdos_data
        ).bandspdosfigure
        y_max = max(self.bands_pdos.data[0].y)
        y_min = min(self.bands_pdos.data[0].y)
        x_max = max(self.bands_pdos.data[1].x)
        self.bands_pdos.update_layout(
            xaxis=dict(title="q-points"),
            yaxis=dict(title="Phonon Bands (THz)", range=[y_min - 0.1, y_max + 0.1]),
            xaxis2=dict(range=[0, x_max + 0.1]),
        )
        self.children = (self.bands_pdos, *self.children)
        self._model.update_thermo_plot(self.thermal_plot)
