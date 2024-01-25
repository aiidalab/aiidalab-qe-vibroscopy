from aiidalab_qe_vibroscopy.app.settings import Setting
from aiidalab_qe_vibroscopy.app.workchain import workchain_and_builder
from aiidalab_qe_vibroscopy.app.result import Result
from aiidalab_qe.common.panel import OutlinePanel

from aiidalab_widgets_base import ComputationalResourcesWidget


class Outline(OutlinePanel):
    title = "Vibrational properties"
    # description = "IR and Raman spectra; you may also select phononic and dielectric properties"


phonopy_code = ComputationalResourcesWidget(
    description="phonopy",
    default_calc_job_plugin="phonopy.phonopy",
)

property = {
    "outline": Outline,
    "code": {"phonopy": phonopy_code},
    "setting": Setting,
    "workchain": workchain_and_builder,
    "result": Result,
}
