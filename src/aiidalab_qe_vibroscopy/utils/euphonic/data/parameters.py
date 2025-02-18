"""Set of parameters for given Euphonic calculation.

We distinguish between parameters for a single crystal and for a powder calculation: the former requires a path in reciprocal space,
while the latter requires a range of q-points.
We have a set of common parameters that are shared between the two types of calculations.
"""

common_parameters = {
    "weighting": "coherent",  # Spectral weighting to plot: DOS, coherent inelastic neutron scattering (default: dos)
    "grid": None,  # FWHM of broadening on q axis in 1/LENGTH_UNIT (no broadening if unspecified). (default: None)
    "grid_spacing": 0.1,  # q-point spacing of Monkhorst-Pack grid. (default: 0.1)
    "energy_unit": "meV",
    "temperature": 0,  # Temperature in K; enable Debye-Waller factor calculation. (Only applicable when --weighting=coherent). (default: None)
    "shape": "gauss",  # The broadening shape (default: gauss)
    "length_unit": "angstrom",
    "q_spacing": 0.01,  # Target distance between q-point samples in 1/LENGTH_UNIT (default: 0.025)
    "energy_broadening": 1,
    "q_broadening": None,  # FWHM of broadening on q axis in 1/LENGTH_UNIT (no broadening if unspecified). (default: None)
    "ebins": 200,  # Number of energy bins (default: 200)
    "e_min": 0,
    "e_max": None,
    "title": None,
    "ylabel": "meV",
    "xlabel": "",
    "save_json": False,
    "no_base_style": False,
    "style": False,
    "vmin": None,
    "vmax": None,
    "save_to": None,
    "asr": None,  # Apply an acoustic-sum-rule (ASR) correction to the data: "realspace" applies the correction to the force constant matrix in real space. "reciprocal" applies the correction to the dynamical matrix at each q-point. (default: None)
    "dipole_parameter": 1.0,  # Set the cutoff in real/reciprocal space for the dipole Ewald sum; higher values use more reciprocal terms. If tuned correctly this can result in performance improvements. See euphonic-optimise-dipole-parameter program for help on choosing a good DIPOLE_PARAMETER. (default: 1.0)
    "use_c": True,
    "n_threads": 1,
}

parameters_single_crystal = {
    **common_parameters,
}

parameters_powder = {
    **common_parameters,
    "q_min": 0,
    "q_max": 1,
    "npts": 150,
    "npts_density": None,
    "pdos": None,
    "e_i": None,
    "sampling": "golden",
    "jitter": True,
    "e_f": None,
    "disable_widgets": True,
}
