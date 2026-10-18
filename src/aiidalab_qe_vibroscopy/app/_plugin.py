from aiidalab_qe.common.panel import PluginOutline

from aiidalab_qe_vibroscopy.app.model import VibroConfigurationSettingsModel
from aiidalab_qe_vibroscopy.app.settings import VibroConfigurationSettingPanel
from aiidalab_qe_vibroscopy.app.code import (
    VibroResourceSettingsModel,
    VibroResourcesSettingsPanel,
)
from aiidalab_qe_vibroscopy.app.result.result import VibroResultsPanel
from aiidalab_qe_vibroscopy.app.result.model import VibroResultsModel

from aiidalab_qe_vibroscopy.app.workchain import workchain_and_builder

from pathlib import Path


class VibroPluginOutline(PluginOutline):
    title = "Vibrational spectroscopy (phonons)"


property = {
    "outline": VibroPluginOutline,
    "configuration": {
        "panel": VibroConfigurationSettingPanel,
        "model": VibroConfigurationSettingsModel,
    },
    "resources": {
        "panel": VibroResourcesSettingsPanel,
        "model": VibroResourceSettingsModel,
    },
    "result": {
        "panel": VibroResultsPanel,
        "model": VibroResultsModel,
    },
    "workchain": workchain_and_builder,
    "guides": {
        "title": "Phonons",
        "path": Path(__file__).resolve().parent / "guides",
    },
}
