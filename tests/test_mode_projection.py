"""Numerical checks for atom mapping and participation-weighted spectra."""

from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from aiida.common.extendeddicts import AttributeDict
from aiida_vibroscopy.data.vibro_mixin import VibrationalMixin

from aiidalab_qe_vibroscopy.app.widgets.ramanmodel import RamanModel
from aiidalab_qe_vibroscopy.utils.mode_projection import (
    atom_participation,
    unitcell_to_primitive,
)
from aiidalab_qe_vibroscopy.utils.phonons.result import atomic_pdos_from_workchain


def make_phonopy(structure, primitive="auto"):
    # Real masses and positions, deliberately dummy species as in AiiDA kinds.
    unitcell = PhonopyAtoms(
        symbols=["He"] * len(structure),
        masses=structure.get_masses(),
        cell=structure.cell,
        positions=structure.positions,
    )
    return Phonopy(unitcell, [2, 1, 1], primitive_matrix=primitive)


def test_mass_convention_and_complex_phase():
    modes = np.array([[[1.0j, 0, 0], [1, 0, 0]]])
    weights = atom_participation(modes, [12, 1], [100])
    np.testing.assert_allclose(weights, [[12 / 13, 1 / 13]])
    np.testing.assert_allclose(
        atom_participation(modes * np.exp(0.7j), [12, 1], [100]), weights
    )


def test_degenerate_subspace_is_rotation_invariant():
    masses = np.array([12.0, 1.0])
    eigenvectors = np.zeros((2, 2, 3))
    eigenvectors[0, 0, 0] = 1
    eigenvectors[1, 1, 0] = 1
    modes = eigenvectors / np.sqrt(masses)[None, :, None]
    rotation = np.array([[0.6, 0.8], [-0.8, 0.6]])
    rotated = np.einsum("ij,jka->ika", rotation, modes)
    original = atom_participation(modes, masses, [100, 100])
    np.testing.assert_allclose(original, [[0.5, 0.5], [0.5, 0.5]])
    np.testing.assert_allclose(
        atom_participation(rotated, masses, [100, 100]), original
    )
    np.testing.assert_allclose(atom_participation(modes, masses, [100, 101]), np.eye(2))


def test_mapping_reduced_cell_and_dummy_species():
    structure = bulk("Al", "fcc", a=4.0, cubic=True)
    phonopy = make_phonopy(structure)
    assert len(phonopy.primitive) == 1
    np.testing.assert_array_equal(
        unitcell_to_primitive(phonopy, structure), [0, 0, 0, 0]
    )
    reordered = structure[[1, 0, 2, 3]]
    with pytest.raises(ValueError, match="numbering"):
        unitcell_to_primitive(phonopy, reordered)


@pytest.mark.parametrize("invalid", [0.0, np.nan, np.inf])
def test_invalid_mode_norm(invalid):
    modes = np.full((1, 2, 3), invalid)
    with pytest.raises(ValueError):
        atom_participation(modes, [12, 1], [100])


def test_atomic_pdos_uses_mapping_and_input_cell_normalization():
    structure = bulk("Al", "fcc", a=4.0, cubic=True)
    phonopy = make_phonopy(structure)
    curve = np.array([0.0, 1.0, 2.0])
    source = SimpleNamespace(get_phonopy_instance=lambda: phonopy)
    creator = SimpleNamespace(
        inputs=AttributeDict(
            {
                "force_constants": source,
                "parameters": SimpleNamespace(get_dict=lambda: {"pdos": "auto"}),
            }
        )
    )
    pdos = SimpleNamespace(creator=creator, get_y=lambda: [("PDOS", curve, "1/THz")])
    workchain = SimpleNamespace(
        inputs=SimpleNamespace(structure=SimpleNamespace(get_ase=lambda: structure)),
        outputs=SimpleNamespace(phonon_pdos=pdos),
    )
    atomic = atomic_pdos_from_workchain(workchain)
    np.testing.assert_allclose(atomic, np.tile(curve, (4, 1)))
    np.testing.assert_allclose(atomic.sum(axis=0), 4 * curve)


class SyntheticProjectedData(VibrationalMixin):
    def __init__(self):
        self.structure = Atoms(
            "CH", positions=[[0, 0, 0], [1, 0, 0]], cell=[4, 5, 6], pbc=True
        )
        self.phonopy = make_phonopy(self.structure, np.eye(3))

    def get_phonopy_instance(self):
        return self.phonopy

    def mode_data(self, nac_direction=None):
        frequencies = np.array([100.0, 400.0])
        weights = np.array([[0.8, 0.2], [0.1, 0.9]])
        if nac_direction is not None:
            frequencies += 10
            weights = weights[:, ::-1]
        displacements = np.zeros((2, 2, 3))
        displacements[:, :, 0] = np.sqrt(weights / self.structure.get_masses())
        return frequencies, displacements, ["A", "B"]

    def run_active_modes(self, selection_rule=None, nac_direction=None):
        return self.mode_data(nac_direction)

    def run_raman_susceptibility_tensors(self, **kwargs):
        frequencies, _, labels = self.mode_data(kwargs.get("nac_direction"))
        return np.array([np.eye(3), np.diag([2.0, 3.0, 4.0])]), frequencies, labels

    def run_polarization_vectors(self, **kwargs):
        frequencies, _, labels = self.mode_data(kwargs.get("nac_direction"))
        return np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]), frequencies, labels


def projected_model(spectrum="Raman", plot_type="powder", separate=False):
    data = SyntheticProjectedData()
    model = RamanModel(
        spectrum_type=spectrum,
        plot_type=plot_type,
        separate_polarizations=separate,
        input_structure=data.structure,
        vibrational_data_uuid="synthetic-loader",
    )
    model.get_raman_data = lambda: data
    model.fetch_data()
    return model


@pytest.mark.parametrize(
    "spectrum, plot_type, separate",
    [
        ("Raman", "powder", False),
        ("Raman", "powder", True),
        ("Raman", "single_crystal", False),
        ("Raman", "plane_average", False),
        ("IR", "powder", False),
        ("IR", "single_crystal", False),
    ],
)
def test_selected_complement_and_all_atoms(spectrum, plot_type, separate):
    model = projected_model(spectrum, plot_type, separate)
    model.selected_atoms = "1"
    model.update_data()
    first = model.projected_intensities.copy()
    first_depolarized = model.projected_depolarized.copy()
    np.testing.assert_allclose(model.projection_weights, [0.8, 0.1])
    model.selected_atoms = "2"
    model.update_data()
    np.testing.assert_allclose(first + model.projected_intensities, model.intensities)
    if separate:
        np.testing.assert_allclose(
            first_depolarized + model.projected_depolarized,
            model.intensities_depolarized,
        )
    model.selected_atoms = "1..2 1"
    model.update_data()
    np.testing.assert_allclose(model.projected_intensities, model.intensities)
    model.selected_atoms = ""
    model.update_data()
    assert model.projected_intensities.size == 0


def test_nac_change_refreshes_mode_participation():
    model = projected_model()
    model.selected_atoms = "1"
    model.update_data()
    np.testing.assert_allclose(model.projection_weights, [0.8, 0.1])
    model.use_nac_direction = True
    model.nac_direction = "1 0 0"
    model.update_data()
    np.testing.assert_allclose(model.projection_weights, [0.2, 0.9])
    np.testing.assert_allclose(model._mode_frequencies, model.raw_frequencies)


def test_selected_only_and_reset_do_not_accumulate_traces():
    import plotly.graph_objects as go

    model = projected_model(separate=True)
    figure = go.Figure()
    for _ in range(3):
        model.selected_atoms = "1"
        model.selected_only = False
        model.update_data()
        model.update_plot(figure)
        assert len(figure.data) == 4
        model.selected_only = True
        model.update_plot(figure)
        assert len(figure.data) == 2
        assert all(trace.line.width == 3.5 for trace in figure.data)
        model.selected_atoms = ""
        model.update_data()
        model.update_plot(figure)
        assert len(figure.data) == 2


def test_selected_export_preserves_definition_and_numbering():
    import base64
    import json

    model = projected_model()
    model.selected_atoms = "2"
    model.update_data()
    captured = {}
    model._download = lambda payload, filename: captured.update(
        data=json.loads(base64.b64decode(payload)), filename=filename
    )
    model.download_data()
    assert captured["data"]["Selected atoms (1-based)"] == [2]
    assert "Mass-weighted" in captured["data"]["Projection definition"]
    np.testing.assert_allclose(captured["data"]["Mode participation"], [0.2, 0.9])
