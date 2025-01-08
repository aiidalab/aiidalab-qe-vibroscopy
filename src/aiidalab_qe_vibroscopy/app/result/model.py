from aiidalab_qe.common.panel import ResultsModel
import traitlets as tl


class VibroResultsModel(ResultsModel):
    title = "Vibroscopy"
    identifier = "vibronic"

    _this_process_label = "VibroWorkChain"

    tab_titles = tl.List([])

    def get_vibro_node(self):
        return self._get_child_outputs()

    def needs_dielectric_tab(self):
        node = self.get_vibro_node()
        if not any(key in node for key in ["iraman", "dielectric", "harmonic"]):
            return False
        return True

    def needs_raman_tab(self):
        node = self.get_vibro_node()
        if not any(key in node for key in ["iraman", "harmonic"]):
            return False
        return True

    # Here we use _fetch_child_process_node() since the function needs the input_structure in inputs
    def needs_phonons_tab(self):
        node = self.get_vibro_node()
        if not any(
            key in node for key in ["phonon_bands", "phonon_thermo", "phonon_pdos"]
        ):
            return False
        return True

    # Here we use _fetch_child_process_node() since the function needs the input_structure in inputs
    def needs_euphonic_tab(self):
        node = self.get_vibro_node()
        if not any(key in node for key in ["phonon_bands"]):
            return False
        return True
