"""Bands results view widgets"""

from __future__ import annotations


import numpy as np


from aiida_vibroscopy.utils.broadenings import multilorentz


def plot_powder(
    frequencies: list[float],
    intensities: list[float],
    broadening: float = 10.0,
    x_range: list[float] | str = "auto",
    broadening_function=multilorentz,
    normalize: bool = True,
):
    frequencies = np.array(frequencies)
    intensities = np.array(intensities)

    if x_range == "auto":
        xi = max(0, frequencies.min() - 200)
        xf = frequencies.max() + 200
        x_range = np.arange(xi, xf, 1.0)

    y_range = broadening_function(x_range, frequencies, intensities, broadening)

    if normalize:
        y_range /= y_range.max()

    return x_range, y_range


def export_iramanworkchain_data(node):
    """
    We have multiple choices: IR, RAMAN.
    """

    if "iraman" in node:
        output_node = node.iraman
    elif "harmonic" in node:
        output_node = node.harmonic
    else:
        # we have raman and ir only if we run IRamanWorkChain or HarmonicWorkChain
        return None

    if "vibrational_data" in output_node:
        # We enable the possibility to provide both spectra.
        # We give as output or string, or the output node.

        spectra_data = {
            "Raman": None,
            "Ir": None,
        }

        vibrational_data = output_node.vibrational_data
        vibro = (
            vibrational_data.numerical_accuracy_4
            if hasattr(vibrational_data, "numerical_accuracy_4")
            else vibrational_data.numerical_accuracy_2
        )

        if "born_charges" in vibro.get_arraynames():
            (
                polarized_intensities,
                frequencies,
                labels,
            ) = vibro.run_powder_ir_intensities()
            total_intensities = polarized_intensities

            # sometimes IR/Raman has not active peaks by symmetry, or due to the fact that 1st order cannot capture them
            if len(total_intensities) == 0:
                spectra_data["Ir"] = (
                    "No IR modes detected."  # explanation added in the main results script of the app.
                )
            else:
                spectra_data["Ir"] = output_node

        if "raman_tensors" in vibro.get_arraynames():
            (
                polarized_intensities,
                depolarized_intensities,
                frequencies,
                labels,
            ) = vibro.run_powder_raman_intensities(frequency_laser=532, temperature=300)
            total_intensities = polarized_intensities + depolarized_intensities

            # sometimes IR/Raman has not active peaks by symmetry, or due to the fact that 1st order cannot capture them
            if len(total_intensities) == 0:
                spectra_data["Raman"] = (
                    "No Raman modes detected."  # explanation added in the main results script of the app.
                )
            else:
                spectra_data["Raman"] = output_node

        return spectra_data
    else:
        return None


def fix_eigenvectors_dimensions(model):
    """Fix eigenvectors dimensions from primitive to supercell if needed."

    This is needed when the vibrational analysis is done on a primitive cell, i.e. the eigenvectors are less than
    the number of atoms in the supercell. This gives issues in the visualization of the modes.
    """

    vibro = (
        model.vibro["iraman"] if "iraman" in model.vibro else model.vibro["harmonic"]
    )
    if vibro is None:
        return None

    if "numerical_accuracy_4" in vibro.get("vibrational_data", {}):
        vibro_data = vibro.vibrational_data.numerical_accuracy_4
    elif "numerical_accuracy_2" in vibro.get("vibrational_data", {}):
        vibro_data = vibro.vibrational_data.numerical_accuracy_2
    else:
        raise ValueError(
            "No vibrational data found in the provided vibro object.", vibro.keys()
        )

    ph_instance = vibro_data.get_phonopy_instance()
    structure = model.input_structure

    shape_primitive = model.eigenvectors.shape
    sc_atomic_numbers = structure.get_atomic_numbers()
    sc_eigenvectors = np.zeros(
        (shape_primitive[0], len(sc_atomic_numbers), shape_primitive[2])
    )

    p2s_map = ph_instance.primitive.p2s_map
    s2p_map = ph_instance.primitive.s2p_map
    if len(p2s_map) == len(s2p_map):
        # nothing to do
        return model.eigenvectors

    for i, t in enumerate(p2s_map.tolist()):
        where = np.where(s2p_map == t)
        print(i, t, where)
        for j in range(shape_primitive[0]):
            sc_eigenvectors[j, where[0], :] = model.eigenvectors[j, i, :]

    return sc_eigenvectors
