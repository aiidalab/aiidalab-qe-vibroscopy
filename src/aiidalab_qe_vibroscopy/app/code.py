from aiidalab_qe.common.code.model import CodeModel, PwCodeModel
from aiidalab_qe.common.panel import (
    PluginResourceSettingsModel,
    PluginResourceSettingsPanel,
)


class VibroResourceSettingsModel(PluginResourceSettingsModel):
    """Resource settings for the vibroscopy calculations."""

    title = "Vibronic resources"
    identifier = "vibronic"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add_models(
            {
                "phonon": PwCodeModel(
                    name="phonon",
                    description="pw.x for phonons",
                    default_calc_job_plugin="quantumespresso.pw",
                ),
                "dielectric": PwCodeModel(
                    name="dielectric",
                    description="pw.x for dielectric",
                    default_calc_job_plugin="quantumespresso.pw",
                ),
                "phonopy": CodeModel(
                    name="phonopy",
                    description="phonopy",
                    default_calc_job_plugin="phonopy.phonopy",
                ),
            }
        )


class VibroResourcesSettingsPanel(
    PluginResourceSettingsPanel[VibroResourceSettingsModel]
):
    """Panel for the resource settings for the vibroscopy calculations."""

    title = "Vibronic"
