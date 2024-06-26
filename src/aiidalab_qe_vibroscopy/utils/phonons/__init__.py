# from aiidalab_qe_vibroscopy.phonons.settings import Setting
# from aiidalab_qe_vibroscopy.phonons.workchain import workchain_and_builder
# from aiidalab_qe_vibroscopy.phonons.result import Result
from aiidalab_qe.common.panel import OutlinePanel
import ipywidgets as ipw


class Outline(OutlinePanel):
    title = "Phonon band structure"
    help = "Harmonic approximation"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.layout = ipw.Layout(width="600px", display="none")


property = {
    "outline": Outline,
    # "setting": Setting,
    # "workchain": workchain_and_builder,
    # "result": Result,
}
