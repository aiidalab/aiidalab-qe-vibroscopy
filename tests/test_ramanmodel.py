"""Test Raman plotting with synthetic tensors, without an AiiDA profile."""

import numpy as np
import pytest
from aiida_vibroscopy.data.vibro_mixin import VibrationalMixin
from aiida_vibroscopy.utils.broadenings import multilorentz

from aiidalab_qe_vibroscopy.app.widgets.ramanmodel import RamanModel


class SyntheticRamanData(VibrationalMixin):
    """Use the backend intensity methods with two known Cartesian tensors."""

    def run_raman_susceptibility_tensors(self, **kwargs):
        return (
            np.array(
                [
                    [[1.0, 0.5, 0.0], [0.5, 2.0, 1.0], [0.0, 1.0, 4.0]],
                    [[4.0, 1.0, 0.5], [1.0, 1.0, 0.0], [0.5, 0.0, 2.0]],
                ]
            ),
            np.array([100.0, 400.0]),
            ["A", "B"],
        )

    def run_polarization_vectors(self, **kwargs):
        return (
            np.array([[1.0, 2.0, 4.0], [4.0, 1.0, 2.0]]),
            np.array([100.0, 400.0]),
            ["A", "B"],
        )


@pytest.mark.parametrize(
    "plane, axes", [("xy", (0, 1)), ("yz", (1, 2)), ("xz", (0, 2))]
)
@pytest.mark.parametrize(
    "frequency_laser, temperature", [(532.0, 300.0), (785.0, 300.0), (532.0, 100.0)]
)
def test_plane_average_matches_cartesian_backend(
    plane, axes, frequency_laser, temperature
):
    """Preserve the backend's Cartesian axes and mode-dependent Raman weights."""
    data = SyntheticRamanData()
    model = RamanModel(
        spectrum_type="Raman",
        plot_type="plane_average",
        plane_type=plane,
        frequency_laser=frequency_laser,
        temperature=temperature,
    )
    model.get_raman_data = lambda: data
    model.update_data()

    # Integrating independent incoming/outgoing angles over [0, 2*pi] gives
    # pi**2 times the sum of the four Cartesian polarization intensities.
    # Use the backend's projection and prefactor, independently of the app.
    expected = np.zeros(2)
    cartesian = np.eye(3)
    for incoming in axes:
        for outgoing in axes:
            intensities, frequencies, _ = data.run_single_crystal_raman_intensities(
                pol_incoming=cartesian[incoming],
                pol_outgoing=cartesian[outgoing],
                frequency_laser=frequency_laser,
                temperature=temperature,
                absolute=True,
            )
            expected += np.pi**2 * intensities

    np.testing.assert_array_equal(model.raw_frequencies, frequencies)
    np.testing.assert_allclose(model.raw_intensities, expected, rtol=1e-10, atol=0)
    # A shared, incorrect scalar prefactor also corrupts relative peak weights.
    np.testing.assert_allclose(
        model.raw_intensities / np.max(model.raw_intensities),
        expected / np.max(expected),
        rtol=1e-10,
        atol=0,
    )


@pytest.mark.parametrize("plot_type", ["powder", "single_crystal"])
@pytest.mark.parametrize("frequency_laser", [532.0, 785.0])
def test_raman_laser_wavelength_matches_backend(plot_type, frequency_laser):
    """Apply the selected laser wavelength in both backend intensity paths."""
    data = SyntheticRamanData()
    model = RamanModel(
        spectrum_type="Raman",
        plot_type=plot_type,
        frequency_laser=frequency_laser,
        temperature=100.0,
    )
    model.get_raman_data = lambda: data
    model.update_data()

    if plot_type == "powder":
        hh, hv, _, _ = data.run_powder_raman_intensities(
            frequency_laser=frequency_laser, temperature=100.0
        )
        expected = hh + hv
    else:
        expected, _, _ = data.run_single_crystal_raman_intensities(
            pol_incoming=[0.0, 0.0, 1.0],
            pol_outgoing=[0.0, 0.0, 1.0],
            frequency_laser=frequency_laser,
            temperature=100.0,
        )

    np.testing.assert_allclose(model.raw_intensities, expected, rtol=1e-10, atol=0)


@pytest.mark.parametrize("spectrum_type", ["Raman", "IR"])
@pytest.mark.parametrize("broadening", [10.0, 40.0])
def test_single_crystal_uses_selected_broadening(spectrum_type, broadening):
    """The single-crystal line shape must follow the selected FWHM."""
    model = RamanModel(
        spectrum_type=spectrum_type,
        plot_type="single_crystal",
        broadening=broadening,
    )
    model.get_raman_data = SyntheticRamanData
    model.update_data()

    expected = multilorentz(
        model.frequencies, model.raw_frequencies, model.raw_intensities, broadening
    )
    expected /= expected.max()
    np.testing.assert_allclose(model.intensities, expected, rtol=1e-10, atol=0)
