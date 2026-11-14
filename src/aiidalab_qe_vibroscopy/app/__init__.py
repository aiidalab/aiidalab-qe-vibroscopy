"""AiiDAlab QE plugin for vibrational spectroscopy (phonons)."""

from aiidalab_qe.plugins import QeAppPlugin


class VibroPlugin(QeAppPlugin):
    @property
    def outline(self) -> str:
        return "Vibrational spectroscopy (phonons)"

    def configuration(self):
        from .model import VibroConfigurationSettingsModel
        from .settings import VibroConfigurationSettingPanel

        return {
            "panel": VibroConfigurationSettingPanel,
            "model": VibroConfigurationSettingsModel,
        }

    def resources(self):
        from .code import VibroResourceSettingsModel, VibroResourcesSettingsPanel

        return {
            "panel": VibroResourcesSettingsPanel,
            "model": VibroResourceSettingsModel,
        }

    def result(self):
        from .result.model import VibroResultsModel
        from .result.result import VibroResultsPanel

        return {
            "panel": VibroResultsPanel,
            "model": VibroResultsModel,
        }

    def workchain(self):
        from .workchain import workchain_and_builder

        return workchain_and_builder

    def guides(self):
        from pathlib import Path

        return {
            "title": "Phonons",
            "path": Path(__file__).resolve().parent / "guides",
        }
