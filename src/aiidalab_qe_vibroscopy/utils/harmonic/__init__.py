# from aiidalab_qe_vibroscopy.harmonic.settings import Setting
# from aiidalab_qe_vibroscopy.harmonic.workchain import workchain_and_builder
# from aiidalab_qe_vibroscopy.harmonic.result import Result
from aiidalab_qe.common.panel import OutlinePanel


class Outline(OutlinePanel):
    title = "Phonon properties"
    # description = "Select to proceed with the calculation of the phononic and dielectric properties"


property = {
    "outline": Outline,
    # "setting": Setting,
    # "workchain": workchain_and_builder,
    # "result": Result,
}
