from html import escape
import ipywidgets as ipw

from aiidalab_widgets_base import LoadingWidget
from aiidalab_qe_vibroscopy.app.widgets.phononmodel import PhononModel
from aiidalab_qe_vibroscopy.app.widgets.atom_selection import AtomSelectionControls

import plotly.graph_objects as go
from aiidalab_qe.common.bands_pdos.bandpdosplotly import BandsPdosPlotly

from aiidalab_qe.common.infobox import InAppGuide


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
        self._model.process_uuid = node if isinstance(node, str) else node.uuid
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
            yaxis=dict(title="THz", range=[y_min - 0.1, y_max + 0.1]),
            xaxis2=dict(range=[0, x_max + 0.1]),
        )
        self.atom_controls = AtomSelectionControls(len(self._model.input_structure))
        self.atom_controls.selected_only.description = "Show selected DOS only"
        self.apply_selection = ipw.Button(
            description="Update atom selection", button_style="primary"
        )
        self.apply_selection.on_click(self._on_atom_selection)
        self._dos_traces = [
            trace for trace in self.bands_pdos.data if trace.xaxis == "x2"
        ]
        self._selected_trace = None
        if self._model.projection_error:
            self.atom_controls.message.value = escape(self._model.projection_error)
        self.children = (
            InAppGuide(identifier="phonons-spectrum-results"),
            self.atom_controls,
            self.apply_selection,
            self.bands_pdos,
            *self.children,
        )
        self._model.update_thermo_plot(self.thermal_plot)

    def _on_atom_selection(self, _=None):
        self._model.selected_atoms = self.atom_controls.selection.value
        self._model.selected_only = self.atom_controls.selected_only.value
        try:
            self._model.update_atom_selection()
        except ValueError as exc:
            self.atom_controls.message.value = (
                f"<div class='alert alert-danger'>{escape(str(exc))}</div>"
            )
            return
        selected = self._model.selected_pdos
        self.atom_controls.message.value = (
            f"{len(self._model.selected_indices)} atoms selected."
            if self._model.selected_indices
            else ""
        )
        with self.bands_pdos.batch_update():
            for trace in self._dos_traces:
                trace.visible = not (
                    self._model.selected_only and self._model.selected_indices
                )
            if selected is None:
                if self._selected_trace is not None:
                    self._selected_trace.visible = False
                return
            if self._selected_trace is None:
                reference = self._dos_traces[0]
                self.bands_pdos.add_trace(
                    go.Scatter(
                        x=selected,
                        y=self._model.pdos_data["dos"][0]["x"],
                        xaxis=reference.xaxis,
                        yaxis=reference.yaxis,
                        name="Selected atoms",
                        line={"width": 3.5, "color": "#d62728"},
                    )
                )
                self._selected_trace = self.bands_pdos.data[-1]
            self._selected_trace.x = selected
            self._selected_trace.visible = True

    def close(self):
        if hasattr(self, "apply_selection"):
            self.apply_selection.on_click(self._on_atom_selection, remove=True)
            self.apply_selection.close()
        if hasattr(self, "atom_controls"):
            self.atom_controls.close()
        super().close()
