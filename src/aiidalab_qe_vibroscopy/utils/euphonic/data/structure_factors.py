from typing import List, Optional


import matplotlib.style
import numpy as np
import copy
import seekpath
from math import ceil

""""
Check double imports!
"""
import euphonic
from euphonic import ureg, QpointFrequencies, ForceConstants
import euphonic.plot
from euphonic.util import get_qpoint_labels
from euphonic.styles import base_style
from euphonic.cli.utils import (
    _bands_from_force_constants,
    _calc_modes_kwargs,
    _compose_style,
    _plot_label_kwargs,
    _get_debye_waller,
    _get_energy_bins,
    _get_q_distance,
    matplotlib_save_or_show,
)

from euphonic.cli.utils import (
    _get_pdos_weighting,
    _arrange_pdos_groups,
)
from euphonic.powder import (
    sample_sphere_dos,
    sample_sphere_pdos,
    sample_sphere_structure_factor,
)
from euphonic.spectra import apply_kinematic_constraints
from euphonic.styles import intensity_widget_style
import euphonic.util


from aiidalab_qe_vibroscopy.utils.euphonic.data.parameters import (
    parameters_single_crystal,
    parameters_powder,
)

# Dummy tqdm function if tqdm progress bars unavailable
try:
    from tqdm.auto import tqdm
except ModuleNotFoundError:
    try:
        from tqdm import tqdm
    except ModuleNotFoundError:

        def tqdm(sequence):
            return sequence


import sys
import os


# Disable
def blockPrint():
    sys.stdout = open(os.devnull, "w")


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


########################
################################ START DESCRIPTION
########################

"""
In this module we have the functions used to obtain the intensity maps (dynamical structure factor)
and the powder maps(from euphonic, using the force constants instances as obtained from phonopy.yaml).
These are then used in the widgets to plot the corresponding quantities.

NOTE: the two main functions here are produce_bands_weigthed_data and produce_powder_data.

PLEASE NOTE: the scattering lengths are tabulated (Euphonic/euphonic/data/sears-1992.json) and are from Sears (1992) Neutron News 3(3) pp26--37.
"""

########################
################################ END DESCRIPTION
########################

########################
################################ START custom lin q path routine
########################


def join_q_paths(coordinates: list, labels: list, delta_q=0.1, G=[0, 0, 0]):
    from euphonic.cli.utils import _get_tick_labels, _get_break_points

    """
    Join list of q points and labels.

    Inputs:

    coordinates (list): list of tuples with the coordinates of path: [(kxi,kyi,kzi),(kxf,kyf,kzf)].
                        Coordinates should be in reciprocal lattice units, as euphonic wants.
    labels (list): list of tuples with the labels: [("Gamma", "M")].
    delta_q (float): q spacing in Angstrom^-1.
    G (list): the modulus of the three reciprocal lattice vectors. Used to convert back and forth into RLU

    Outputs:

    final_path [list]: list of q points with labels.
    refined_labels, split_args: list needed for the produce_bands_weigthed_data method.
    """
    list_of_paths = []
    Nq_tot = 0

    new_labels_index = []  # here we store the index to then set the labels list as in seekpath, to be refined in the produce_curated_data routine.
    for path in coordinates:
        kxi, kyi, kzi = path[0]  # starting point
        kxf, kyf, kzf = path[1]  # end point

        Nq = int(
            np.linalg.norm(
                np.array([(kxf - kxi) * G[0], (kyf - kyi) * G[1], (kzf - kzi) * G[2]])
            )
            / delta_q
        )

        kpath = np.array(
            [
                np.linspace(kxi, kxf, Nq),
                np.linspace(kyi, kyf, Nq),
                np.linspace(kzi, kzf, Nq),
            ]
        )
        list_of_paths.append(kpath.T.reshape(Nq, 3))

        Nq_tot += Nq
        print(Nq_tot, Nq)
        new_labels_index.append(Nq_tot - Nq)
        new_labels_index.append(Nq_tot - 1)

    print(new_labels_index, labels)
    final_path = list_of_paths[0]
    for linear_path in list_of_paths[1:]:
        final_path = np.vstack([final_path, linear_path])

    # final_path = final_path.reshape(len(list_of_paths)*Nq,3)

    # like in _get_tick_labels(bandpath: dict) -> List[Tuple[int, str]] of Euphonic
    new_labels = [""] * len(final_path)
    for ind, label in list(zip(new_labels_index, labels)):
        new_labels[ind] = label
        # print(ind,label)

    bandpath = {"explicit_kpoints_labels": new_labels}
    refined_labels = _get_tick_labels(bandpath={"explicit_kpoints_labels": new_labels})
    split_args = {"indices": _get_break_points(bandpath)}

    return final_path, refined_labels, split_args


########################
################################ END custom q path routine
########################


########################
################################ START INTENSITY PLOT GENERATOR
########################


def produce_bands_weigthed_data(
    params: Optional[List[str]] = None,
    fc: ForceConstants = None,
    linear_path=None,
    plot=False,
) -> None:
    blockPrint()
    """
    This is essentially an adapted version of the function ("main") implemented in Euphonic
    for the cli plotting of the weighted bands. For weighted bands I means or the dynamical structure
    factor or the DOS-weighted one. This can be triggered in inputs,
    and will call the calculate_sqw_map or the calculate_dos_map functions, respectively.

    linear_path is for custom path in the app:
    linear_path = {
        'coordinates':[[(0,0,0),(0.5,0.5,0.5)],[(0.5,0.5,0.5),(1,1,1)]],
        'labels' : ["$\Gamma$","X","X","(1,1,1)"],
        'delta_q':0.1, # A^-1
    }
    """
    # args = get_args(get_parser(), params)
    if not params:
        args = AttrDict(copy.deepcopy(parameters_single_crystal))
    else:
        args = AttrDict(params)

    # redundancy with args...
    calc_modes_kwargs = _calc_modes_kwargs(args)

    frequencies_only = args.weighting != "coherent"
    """data = load_data_from_file(args.filename, verbose=True,
                               frequencies_only=frequencies_only)"""
    data = fc

    if not frequencies_only and type(data) is QpointFrequencies:
        raise TypeError(
            "Eigenvectors are required to use " '"--weighting coherent" option'
        )
    if args.weighting.lower() == "coherent" and args.temperature is not None:
        if not isinstance(data, ForceConstants):
            raise TypeError(
                "Force constants data is required to generate "
                'the Debye-Waller factor. Leave "--temperature" '
                "unset if plotting precalculated phonon modes."
            )

    q_spacing = _get_q_distance(args.length_unit, args.q_spacing)
    recip_length_unit = q_spacing.units

    if isinstance(data, ForceConstants):
        # print("Getting band path...")
        # HERE we add the custom path generation:
        if linear_path:
            # 1. get the rl_norm list for conversion delta_q ==> Nq in the join_q_paths
            structure = fc.crystal.to_spglib_cell()
            bandpath = seekpath.get_explicit_k_path(structure)

            rl = bandpath["reciprocal_primitive_lattice"]
            rl_norm = []
            for G in range(3):
                rl_norm.append(np.linalg.norm(np.array(rl[G])))

            # 2. compute the path via delta_q
            (qpts, x_tick_labels, split_args) = join_q_paths(
                coordinates=linear_path["coordinates"],
                labels=linear_path["labels"],
                delta_q=linear_path["delta_q"],
                G=rl_norm,
            )

            # 3. compute the corresponding phonons
            modes = fc.calculate_qpoint_phonon_modes(
                qpts,
                reduce_qpts=False,
            )

        else:
            # Use seekpath.
            (modes, x_tick_labels, split_args) = _bands_from_force_constants(
                data,
                q_distance=q_spacing,
                # insert_gamma=False,
                insert_gamma=True,
                frequencies_only=frequencies_only,
                **calc_modes_kwargs,
            )
    else:
        modes = data
        x_tick_labels = get_qpoint_labels(
            modes.qpts, cell=modes.crystal.to_spglib_cell()
        )

    # duplication from euphonic/cli/utils.py
    if args.e_min is None:
        # Subtract small amount from min frequency - otherwise due to unit
        # conversions binning of this frequency can vary with different
        # architectures/lib versions, making it difficult to test
        emin_room = 1e-5 * ureg("meV").to(modes.frequencies.units).magnitude
        args.e_min = min(np.min(modes.frequencies.magnitude - emin_room), 0.0)
    if args.e_max is None:
        args.e_max = np.max(modes.frequencies.magnitude) * 1.05
    if args.e_min >= args.e_max:
        raise ValueError(
            f"Maximum energy ({args.e_max}) should be greater than minimum ({args.e_min}). "
        )

    modes.frequencies_unit = args.energy_unit
    ebins = _get_energy_bins(modes, args.ebins + 1, emin=args.e_min, emax=args.e_max)

    # print("Computing intensities and generating 2D maps")

    if args.weighting.lower() == "coherent":
        if args.temperature is not None:
            temperature = args.temperature * ureg("K")
            dw = _get_debye_waller(
                temperature,
                data,
                grid=args.grid,
                grid_spacing=(args.grid_spacing * recip_length_unit),
                **calc_modes_kwargs,
            )
        else:
            dw = None

        spectrum = modes.calculate_structure_factor(dw=dw).calculate_sqw_map(ebins)

    elif args.weighting.lower() == "dos":
        spectrum = modes.calculate_dos_map(ebins)

    if args.q_broadening or args.energy_broadening:
        spectrum = spectrum.broaden(
            x_width=(
                args.q_broadening * recip_length_unit if args.q_broadening else None
            ),
            y_width=(
                args.energy_broadening * ebins.units if args.energy_broadening else None
            ),
            shape=args.shape,
            method="convolve",
        )

    # print("Plotting figure")
    plot_label_kwargs = _plot_label_kwargs(
        args, default_ylabel=f"Energy / {spectrum.y_data.units:~P}"
    )

    if x_tick_labels:
        spectrum.x_tick_labels = x_tick_labels

    spectra = spectrum  # .split(**split_args)  # type: List[Spectrum2D]
    # if len(spectra) > 1:
    #    pass  # print(f"Found {len(spectra)} regions in q-point path")

    if args.save_json:
        spectrum.to_json_file(args.save_json)
    style = _compose_style(user_args=args, base=[base_style])
    if plot:
        with matplotlib.style.context(style):
            fig = euphonic.plot.plot_2d(  # noqa F841
                spectra,  # noqa F841
                vmin=args.vmin,  # noqa F841
                vmax=args.vmax,  # noqa F841
                **plot_label_kwargs,  # noqa F841
            )  # noqa F841
            matplotlib_save_or_show(save_filename=args.save_to)

    enablePrint()

    return spectra, copy.deepcopy(params)


########################
################################ START POWDER
########################


# parameters_powder = AttrDict(par_dict_powder)


def produce_powder_data(
    params: Optional[List[str]] = None,
    fc: ForceConstants = None,
    plot=False,
    linear_path=None,
) -> None:
    blockPrint()
    """Read the description of the produce_bands_weigthed_data function for more details.
    """

    if not params:
        args = AttrDict(copy.deepcopy(parameters_powder))
    else:
        args = AttrDict(params)

    # redundancy with args
    calc_modes_kwargs = _calc_modes_kwargs(args)

    # Make sure we get an error if accessing NPTS inappropriately
    if args.npts_density is not None:
        args.npts = None

    """fc = load_data_from_file(args.filename, verbose=True)"""

    if not isinstance(fc, ForceConstants):
        raise TypeError(
            "Force constants are required to use the " "euphonic-powder-map tool"
        )
    if args.pdos is not None and args.weighting == "coherent":
        raise ValueError(
            '"--pdos" is only compatible with ' '"--weighting" options that include dos'
        )
    # print("Setting up dimensions...")

    q_min = _get_q_distance(args.length_unit, args.q_min)
    q_max = _get_q_distance(args.length_unit, args.q_max)
    recip_length_unit = q_min.units

    n_q_bins = ceil((args.q_max - args.q_min) / args.q_spacing)
    q_bin_edges = (
        np.linspace(q_min.magnitude, q_max.magnitude, n_q_bins + 1, endpoint=True)
        * recip_length_unit
    )
    q_bin_centers = (q_bin_edges[:-1] + q_bin_edges[1:]) / 2

    # Use X-point modes to estimate frequency range, set up energy bins
    # (Not Gamma in case there are only 3 branches; value would be zero!)
    modes = fc.calculate_qpoint_frequencies(
        np.array([[0.0, 0.0, 0.5]]), **calc_modes_kwargs
    )
    modes.frequencies_unit = args.energy_unit

    if args.e_i is not None and args.e_max is None:
        emax = args.e_i
    else:
        emax = args.e_max

    energy_bins = _get_energy_bins(
        modes, args.ebins + 1, emin=args.e_min, emax=emax, headroom=1.2
    )  # Generous headroom as we only checked one q-point

    if args.weighting in ("coherent",):
        # Compute Debye-Waller factor once for re-use at each mod(q)
        # (If temperature is not set, this will be None.)
        if args.temperature is not None:
            temperature = args.temperature * ureg("K")
            dw = _get_debye_waller(
                temperature,
                fc,
                grid=args.grid,
                grid_spacing=(args.grid_spacing * recip_length_unit),
                **calc_modes_kwargs,
            )
        else:
            temperature = None
            dw = None

    # print(f"Sampling {n_q_bins} |q| shells between {q_min:~P} and {q_max:~P}")

    z_data = np.empty((n_q_bins, len(energy_bins) - 1))

    # for q_index in tqdm(range(n_q_bins)):
    for q_index in range(n_q_bins):
        q = q_bin_centers[q_index]

        if args.npts_density is not None:
            npts = ceil(args.npts_density * (q / recip_length_unit) ** 2)
            npts = max(args.npts_min, min(args.npts_max, npts))
        else:
            npts = args.npts

        if args.weighting == "dos" and args.pdos is None:
            spectrum_1d = sample_sphere_dos(
                fc,
                q,
                npts=npts,
                sampling=args.sampling,
                jitter=args.jitter,
                energy_bins=energy_bins,
                **calc_modes_kwargs,
            )
        elif "dos" in args.weighting:
            spectrum_1d_col = sample_sphere_pdos(
                fc,
                q,
                npts=npts,
                sampling=args.sampling,
                jitter=args.jitter,
                energy_bins=energy_bins,
                weighting=_get_pdos_weighting(args.weighting),
                **calc_modes_kwargs,
            )
            spectrum_1d = _arrange_pdos_groups(spectrum_1d_col, args.pdos)
        elif args.weighting == "coherent":
            spectrum_1d = sample_sphere_structure_factor(
                fc,
                q,
                dw=dw,
                temperature=temperature,
                sampling=args.sampling,
                jitter=args.jitter,
                npts=npts,
                energy_bins=energy_bins,
                **calc_modes_kwargs,
            )

        z_data[q_index, :] = spectrum_1d.y_data.magnitude

    # print(f"Final npts: {npts}")

    spectrum = euphonic.Spectrum2D(
        q_bin_edges, energy_bins, z_data * spectrum_1d.y_data.units
    )

    if args.q_broadening or args.energy_broadening:
        spectrum = spectrum.broaden(
            x_width=(
                args.q_broadening * recip_length_unit if args.q_broadening else None
            ),
            y_width=(
                args.energy_broadening * energy_bins.units
                if args.energy_broadening
                else None
            ),
            shape=args.shape,
        )

    if not (args.e_i is None and args.e_f is None):
        # print("Applying kinematic constraints")
        energy_unit = args.energy_unit
        e_i = args.e_i * ureg(energy_unit) if (args.e_i is not None) else None
        e_f = args.e_f * ureg(energy_unit) if (args.e_f is not None) else None
        spectrum = apply_kinematic_constraints(
            spectrum, e_i=e_i, e_f=e_f, angle_range=args.angle_range
        )

    # print(f"Plotting figure: max intensity "
    # f"{np.nanmax(spectrum.z_data.magnitude) * spectrum.z_data.units:~P}")
    plot_label_kwargs = _plot_label_kwargs(
        args,
        default_xlabel=f"|q| / {q_min.units:~P}",
        default_ylabel=f"Energy / {spectrum.y_data.units:~P}",
    )

    if args.save_json:
        spectrum.to_json_file(args.save_json)
    if args.disable_widgets:
        base = [base_style]
    else:
        base = [base_style, intensity_widget_style]
    style = _compose_style(user_args=args, base=base)
    if plot:
        with matplotlib.style.context(style):
            fig = euphonic.plot.plot_2d(
                spectrum, vmin=args.vmin, vmax=args.vmax, **plot_label_kwargs
            )

            if args.disable_widgets is False:
                # TextBox only available from mpl 2.1.0
                try:
                    from matplotlib.widgets import TextBox
                except ImportError:
                    args.disable_widgets = True

            if args.disable_widgets is False:
                min_label = f"Min Intensity ({spectrum.z_data.units:~P})"
                max_label = f"Max Intensity ({spectrum.z_data.units:~P})"
                boxw = 0.15
                boxh = 0.05
                x0 = 0.1 + len(min_label) * 0.01
                y0 = 0.025
                axmin = fig.add_axes([x0, y0, boxw, boxh])
                axmax = fig.add_axes([x0, y0 + 0.075, boxw, boxh])
                image = fig.get_axes()[0].images[0]
                cmin, cmax = image.get_clim()
                pad = 0.05
                fmt_str = ".2e" if cmax < 0.1 else ".2f"
                minbox = TextBox(
                    axmin, min_label, initial=f"{cmin:{fmt_str}}", label_pad=pad
                )
                maxbox = TextBox(
                    axmax, max_label, initial=f"{cmax:{fmt_str}}", label_pad=pad
                )

                def update_min(min_val):
                    image.set_clim(vmin=float(min_val))
                    fig.canvas.draw()

                def update_max(max_val):
                    image.set_clim(vmax=float(max_val))
                    fig.canvas.draw()

                minbox.on_submit(update_min)
                maxbox.on_submit(update_max)

    enablePrint()
    return spectrum, copy.deepcopy(params)
    matplotlib_save_or_show(save_filename=args.save_to)


def produce_Q_section_modes(
    fc,
    h,
    k,
    Q0=np.array([0, 0, 0]),
    n_h=100,
    n_k=100,
    h_extension=1,
    k_extension=1,
    temperature=0,
):
    from euphonic import ureg

    # see get_Q_section
    # h: array vector
    # k: array vector
    # Q0: "point" in Q-space used to build the portion of plane, using also the two vectors h and k.
    # n_h, n_k: number of points along the two directions. or better, the two vectors.

    def get_Q_section(h, k, Q0, n_h, n_k, h_extension, k_extension):
        # every point in the space is Q=Q0+dv1*h+dv2*k, which

        q_list = []
        h_list = []
        k_list = []

        for dv1 in np.linspace(-h_extension, h_extension, n_h):
            for dv2 in np.linspace(-k_extension, k_extension, n_k):
                p = Q0 + dv1 * h + dv2 * k
                q_list.append(p)  # Q list
                h_list.append(dv1)  # *h[0])
                k_list.append(dv2)  # *k[1])

        return np.array(q_list), np.array(h_list), np.array(k_list)

    q_array, h_array, k_array = get_Q_section(
        h, k, Q0, n_h + 1, n_k + 1, h_extension, k_extension
    )

    modes = fc.calculate_qpoint_phonon_modes(qpts=q_array, reduce_qpts=False)

    if temperature > 0:
        blockPrint()
        dw = _get_debye_waller(
            temperature * ureg("K"),
            fc,
            # grid_spacing=(args.grid_spacing * recip_length_unit),
            # **calc_modes_kwargs,
        )
        enablePrint()
    else:
        dw = None

    labels = {
        "q": f"Q0={[np.round(i,3) for i in Q0]}",
        "h": f"h={[np.round(i,3) for i in h]}",
        "k": f"k={[np.round(i,3) for i in k]}",
    }

    return modes, q_array, h_array, k_array, labels, dw


def produce_Q_section_spectrum(
    modes,
    q_array,
    h_array,
    k_array,
    ecenter,
    deltaE=0.5,
    bins=10,
    spectrum_type="coherent",
    dw=None,
    labels=None,
):
    from aiidalab_qe_vibroscopy.utils.euphonic.data.structure_factors import (
        blockPrint,
        enablePrint,
    )

    # bins = 10 # hard coded beacuse here it does not change anything.
    ebins = _get_energy_bins(
        modes, bins + 1, emin=ecenter - deltaE, emax=ecenter + deltaE
    )

    blockPrint()
    if (
        spectrum_type == "coherent"
    ):  # Temperature?? For now let's drop it otherwise it is complicated.
        spectrum = modes.calculate_structure_factor(dw=dw).calculate_sqw_map(ebins)
    elif spectrum_type == "dos":
        spectrum = modes.calculate_dos_map(ebins)

    mu = ecenter
    sigma = (deltaE) / 2

    # Gaussian weights.
    weights = np.exp(-((spectrum.y_data.magnitude - mu) ** 2) / 2 * sigma**2) / np.sqrt(
        2 * np.pi * sigma**2
    )
    av_spec = np.average(spectrum.z_data.magnitude, axis=1, weights=weights[:-1])
    enablePrint()

    return av_spec, q_array, h_array, k_array, labels


def generated_curated_data(spectra):
    # here we concatenate the bands groups and create the ticks and labels.

    ticks_positions = []
    ticks_labels = []

    final_xspectra = spectra.x_data.magnitude
    final_zspectra = spectra.z_data.magnitude

    for k in spectra.x_tick_labels:
        ticks_positions.append(k[0])
        # ticks_labels.append("Gamma") if k[1] == '$\\Gamma$' else ticks_labels.append(k[1])
        ticks_labels.append(k[1])

        if len(ticks_positions) > 1:
            if (
                ticks_positions[-1] == ticks_positions[-2] + 1
                and ticks_labels[-1] != ticks_labels[-2]
            ):
                ticks_labels[-2] = ticks_labels[-2] + "|" + ticks_labels[-1]
                ticks_positions.pop()
                ticks_labels.pop()

            if ticks_labels[-1] == ticks_labels[-2]:
                ticks_positions.pop()
                ticks_labels.pop()

    return final_xspectra, final_zspectra, ticks_positions, ticks_labels
