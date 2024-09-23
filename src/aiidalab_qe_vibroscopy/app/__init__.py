from aiidalab_qe_vibroscopy.app.settings import Setting
from aiidalab_qe_vibroscopy.app.workchain import workchain_and_builder
from aiidalab_qe_vibroscopy.app.result import Result
from aiidalab_qe.common.panel import OutlinePanel

from aiidalab_qe.common.widgets import (
    QEAppComputationalResourcesWidget,
    PwCodeResourceSetupWidget,
)


class Outline(OutlinePanel):
    title = "Vibrational properties"
    # description = "IR and Raman spectra; you may also select phononic and dielectric properties"


PhononWorkChainPwCode = PwCodeResourceSetupWidget(
    description="pw.x for phonons",  # code for the PhononWorkChain workflow",
    default_calc_job_plugin="quantumespresso.pw",
)

# The finite electric field does not support npools (does not work with npools>1), so we just set it as QEAppComputationalResourcesWidget
DielectricWorkChainPwCode = QEAppComputationalResourcesWidget(
    description="pw.x for dielectric",  # code for the DielectricWorChain workflow",
    default_calc_job_plugin="quantumespresso.pw",
)

PhonopyCalculationCode = QEAppComputationalResourcesWidget(
    description="phonopy",  # code for the PhonopyCalculation calcjob",
    default_calc_job_plugin="phonopy.phonopy",
)

property = {
    "outline": Outline,
    "code": {
        "phonon": PhononWorkChainPwCode,
        "dielectric": DielectricWorkChainPwCode,
        "phonopy": PhonopyCalculationCode,
    },
    "setting": Setting,
    "workchain": workchain_and_builder,
    "result": Result,
}
