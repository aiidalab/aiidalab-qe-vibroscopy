"""Vibronic results view widgets"""

from aiidalab_qe_vibroscopy.app.result.model import VibroResultsModel
from aiidalab_qe.common.panel import ResultsPanel

from aiidalab_qe_vibroscopy.app.widgets.dielectricwidget import DielectricWidget
from aiidalab_qe_vibroscopy.app.widgets.dielectricmodel import DielectricModel
import ipywidgets as ipw

from aiidalab_qe_vibroscopy.app.widgets.ir_ramanwidget import IRRamanWidget
from aiidalab_qe_vibroscopy.app.widgets.ir_ramanmodel import IRRamanModel

from aiidalab_qe_vibroscopy.app.widgets.phononwidget import PhononWidget
from aiidalab_qe_vibroscopy.app.widgets.phononmodel import PhononModel

from aiidalab_qe_vibroscopy.app.widgets.euphonicwidget import EuphonicWidget
from aiidalab_qe_vibroscopy.app.widgets.euphonicmodel import (
    EuphonicResultsModel as EuphonicModel,
)


class VibroResultsPanel(ResultsPanel[VibroResultsModel]):
    title = "Vibronic"
    identifier = "vibronic"
    workchain_labels = ["vibro"]

    def render(self):
        if self.rendered:
            return

        self.tabs = ipw.Tab(
            layout=ipw.Layout(min_height="250px"),
            selected_index=None,
        )
        self.tabs.observe(
            self._on_tab_change,
            "selected_index",
        )

        tab_data = []
        vibro_node = self._model.get_vibro_node()

        needs_phonons_tab = self._model.needs_phonons_tab()
        if needs_phonons_tab:
            vibroscopy_node = self._model.fetch_child_process_node()
            phonon_model = PhononModel()
            phonon_widget = PhononWidget(
                model=phonon_model,
                node=vibroscopy_node,
            )
            tab_data.append(("Phonons", phonon_widget))

        needs_raman_tab = self._model.needs_raman_tab()
        if needs_raman_tab:
            vibroscopy_node = self._model.fetch_child_process_node()
            input_structure = vibroscopy_node.inputs.structure.get_ase()
            irraman_model = IRRamanModel()
            irraman_widget = IRRamanWidget(
                model=irraman_model,
                node=vibro_node,
                input_structure=input_structure,
            )

            tab_data.append(("Raman/IR spectra", irraman_widget))

        needs_dielectric_tab = self._model.needs_dielectric_tab()
        if needs_dielectric_tab:
            dielectric_model = DielectricModel()
            dielectric_widget = DielectricWidget(
                model=dielectric_model,
                node=vibro_node,
            )
            tab_data.append(("Dielectric Properties", dielectric_widget))

        needs_euphonic_tab = self._model.needs_euphonic_tab()
        if needs_euphonic_tab:
            euphonic_model = EuphonicModel()
            euphonic_widget = EuphonicWidget(
                model=euphonic_model,
                node=vibro_node,
            )
            tab_data.append(("Neutron scattering", euphonic_widget))

        # Assign children and titles dynamically
        self.tabs.children = [content for _, content in tab_data]

        for index, (title, _) in enumerate(tab_data):
            self.tabs.set_title(index, title)

        self.children = [self.tabs]
        self.rendered = True
        self.tabs.selected_index = 0

    def _on_tab_change(self, change):
        if (tab_index := change["new"]) is None:
            return
        self.tabs.children[tab_index].render()  # type: ignore
