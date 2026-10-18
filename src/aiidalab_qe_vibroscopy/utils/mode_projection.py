"""Atom mapping and mass-weighted participation for vibrational viewers."""

import numpy as np
from aiida_vibroscopy.common import UNITS


def unitcell_to_primitive(phonopy, structure):
    """Map displayed unit-cell indices to primitive eigenvector indices.

    AiiDA kinds may be encoded as dummy chemical species in Phonopy. Check
    positions and the actual masses, rather than interpreting those symbols.
    """
    unitcell = phonopy.unitcell
    if len(structure) != len(unitcell):
        raise ValueError("The displayed atoms do not match the vibrational unit cell.")
    if not np.allclose(structure.cell, unitcell.cell, atol=1e-7, rtol=0):
        raise ValueError("The displayed cell does not match the vibrational unit cell.")
    delta = np.linalg.solve(
        np.asarray(unitcell.cell).T,
        (structure.positions - unitcell.positions).T,
    ).T
    periodic = np.asarray(structure.pbc)
    delta[:, periodic] -= np.rint(delta[:, periodic])
    if not np.allclose(delta @ unitcell.cell, 0, atol=1e-6, rtol=0):
        raise ValueError("Cannot match displayed atom numbering to the stored modes.")
    if not np.allclose(structure.get_masses(), unitcell.masses, atol=1e-6, rtol=0):
        raise ValueError("Displayed atomic masses do not match the stored modes.")
    primitive = phonopy.primitive
    mapping = np.array(
        [
            primitive.p2p_map[primitive.s2p_map[index]]
            for index in phonopy.supercell.u2s_map
        ],
        dtype=int,
    )
    if len(mapping) != len(structure) or set(mapping) != set(range(len(primitive))):
        raise ValueError("The primitive-to-unit-cell atom mapping is incomplete.")
    return mapping


def atom_participation(displacements, masses, frequencies):
    """Return mass-weighted per-atom fractions, one row per mode.

    The backend returns u=e/sqrt(M). Thus M*abs(u)**2 recovers the
    dynamical-matrix eigenvector weights used for phonon PDOS.
    Average weights within degenerate subspaces so the displayed optical
    participation does not depend on arbitrary eigenvector rotations.
    These are mode-character weights, not atomic optical intensities.
    """
    displacements = np.asarray(displacements)
    masses = np.asarray(masses, dtype=float)
    frequencies = np.asarray(frequencies)
    if displacements.shape != (len(frequencies), len(masses), 3):
        raise ValueError("Mode vectors, frequencies and atoms have inconsistent sizes.")
    if (
        not np.all(np.isfinite(displacements))
        or not np.all(np.isfinite(masses))
        or not np.all(masses > 0)
        or not np.all(np.isfinite(frequencies))
    ):
        raise ValueError(
            "Modes, frequencies and masses must be finite, with positive masses."
        )
    weights = np.sum(np.abs(displacements) ** 2, axis=2) * masses
    norms = weights.sum(axis=1, keepdims=True)
    if np.any(norms <= 0) or not np.all(np.isfinite(norms)):
        raise ValueError("Cannot project a mode with zero or non-finite norm.")
    weights /= norms
    # Same frequency tolerance as the backend's Gamma eigenvector gauge fixing.
    tolerance = 1e-5 * UNITS.thz_to_cm
    order = np.argsort(frequencies)
    start = 0
    while start < len(order):
        stop = start + 1
        while (
            stop < len(order)
            and abs(frequencies[order[stop]] - frequencies[order[start]]) < tolerance
        ):
            stop += 1
        group = order[start:stop]
        weights[group] = weights[group].mean(axis=0)
        start = stop
    return weights
