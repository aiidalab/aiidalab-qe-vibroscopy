"""Tests for converting Phonopy data for Euphonic."""

from unittest.mock import Mock

from aiidalab_qe_vibroscopy.utils.euphonic.data import phonopy_interface


class AttributeDict(dict):
    """Provide attribute access to the mock calculation inputs."""

    __getattr__ = dict.__getitem__


def test_force_constants_use_current_phonopy_yaml_api(monkeypatch):
    """Serialize a Phonopy instance using its supported YAML conversion method."""
    phonopy = Mock(
        spec=["force_constants", "produce_force_constants", "to_phonopy_yaml"]
    )
    phonopy.force_constants = object()
    phonopy.to_phonopy_yaml.return_value = "phonopy yaml"

    phonopy_data = Mock()
    phonopy_data.get_phonopy_instance.return_value = phonopy
    phonopy_data.get_cells_mappings.return_value = {"primitive": {"p2s_map": [0]}}
    calculation = Mock(
        inputs=AttributeDict(phonopy_data=phonopy_data),
    )
    expected_force_constants = object()

    monkeypatch.setattr(phonopy_interface, "blockPrint", lambda: None)
    monkeypatch.setattr(phonopy_interface, "enablePrint", lambda: None)
    monkeypatch.setattr(phonopy_interface, "write_force_constants_to_hdf5", Mock())
    monkeypatch.setattr(
        phonopy_interface,
        "generate_force_constant_instance_temporary_fix",
        Mock(return_value=expected_force_constants),
    )

    result = phonopy_interface.generate_force_constant_from_phonopy(calculation)

    phonopy.to_phonopy_yaml.assert_called_once_with()
    assert result is expected_force_constants


def test_force_constants_use_legacy_phonopy_yaml_api(monkeypatch):
    """Keep compatibility with Phonopy releases exposing set_phonon_info."""
    phonopy = Mock(spec=["force_constants", "produce_force_constants"])
    phonopy.force_constants = object()

    phonopy_data = Mock()
    phonopy_data.get_phonopy_instance.return_value = phonopy
    phonopy_data.get_cells_mappings.return_value = {"primitive": {"p2s_map": [0]}}
    calculation = Mock(
        inputs=AttributeDict(phonopy_data=phonopy_data),
    )
    expected_force_constants = object()
    yaml_instances = []

    class LegacyPhonopyYaml:
        def __init__(self):
            yaml_instances.append(self)

        def set_phonon_info(self, phonopy_instance):
            self.phonopy = phonopy_instance

        def __str__(self):
            return "phonopy yaml"

    monkeypatch.setattr(phonopy_interface, "blockPrint", lambda: None)
    monkeypatch.setattr(phonopy_interface, "enablePrint", lambda: None)
    monkeypatch.setattr(phonopy_interface, "write_force_constants_to_hdf5", Mock())
    monkeypatch.setattr(
        phonopy_interface,
        "generate_force_constant_instance_temporary_fix",
        Mock(return_value=expected_force_constants),
    )
    monkeypatch.setattr("phonopy.interface.phonopy_yaml.PhonopyYaml", LegacyPhonopyYaml)

    result = phonopy_interface.generate_force_constant_from_phonopy(calculation)

    assert yaml_instances[0].phonopy is phonopy
    assert result is expected_force_constants
