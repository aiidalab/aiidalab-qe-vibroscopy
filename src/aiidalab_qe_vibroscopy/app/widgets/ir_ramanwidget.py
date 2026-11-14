from aiidalab_qe_vibroscopy.app.widgets.ir_ramanmodel import IRRamanModel
from aiidalab_qe_vibroscopy.app.widgets.ramanwidget import RamanWidget
from aiidalab_qe_vibroscopy.app.widgets.ramanmodel import RamanModel
import ipywidgets as ipw

from aiidalab_qe.common.infobox import InAppGuide


class IRRamanWidget(ipw.VBox):
    def __init__(self, model: IRRamanModel, node: None, input_structure, **kwargs):
        super().__init__(
            children=[ipw.HTML("Loading Raman data...")],
            **kwargs,
        )
        self._model = model
        self._model.vibro = node
        self._model.input_structure = input_structure
        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.children = [InAppGuide(identifier="Raman-spectrum-results")]

        self.rendered = True
        self._model.fetch_data()
        self._needs_raman_widget()
        self._needs_ir_widget()
        self.render_widgets()

    def _needs_raman_widget(self):
        if self._model.needs_raman_tab:
            self.raman_model = RamanModel()
            self.raman_widget = RamanWidget(
                model=self.raman_model,
                node=self._model.vibro,
                input_structure=self._model.input_structure,
                spectrum_type="Raman",
            )
            self.children = (*self.children, self.raman_widget)

    def _needs_ir_widget(self):
        if self._model.needs_ir_tab:
            self.ir_model = RamanModel()
            self.ir_widget = RamanWidget(
                model=self.ir_model,
                node=self._model.vibro,
                input_structure=self._model.input_structure,
                spectrum_type="IR",
            )
            self.children = (*self.children, self.ir_widget)
        else:
            self.children = (*self.children, ipw.HTML("No IR modes detected."))

    def render_widgets(self):
        if self._model.needs_raman_tab:
            self.raman_widget.render()
        if self._model.needs_ir_tab:
            self.ir_widget.render()
