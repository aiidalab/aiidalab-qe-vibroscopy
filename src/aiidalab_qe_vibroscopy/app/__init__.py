"""Load the QE plugin only when requested, keeping spectrum widgets reusable."""

from importlib import import_module

__all__ = [
    "PluginOutline",
    "VibroConfigurationSettingsModel",
    "VibroConfigurationSettingPanel",
    "VibroResourceSettingsModel",
    "VibroResourcesSettingsPanel",
    "VibroResultsPanel",
    "VibroResultsModel",
    "workchain_and_builder",
    "Path",
    "VibroPluginOutline",
    "property",
]


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module("._plugin", __name__), name)
    globals()[name] = value
    return value
