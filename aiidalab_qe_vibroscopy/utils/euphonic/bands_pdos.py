import numpy as np

import euphonic as eu
import euphonic.util as util
import euphonic.plot as plt

from aiidalab_qe.common.bandpdoswidget import cmap

import json
from monty.json import jsanitize

import seekpath

from aiida import orm

########################
################################ START DESCRIPTION
########################

"""
In this module we have the functions used to obtain the phonon bands and pdos
(from euphonic, using the force constants instances as obtained from phonopy.yaml).
These are then used in the widgets to plot the corresponding quantities.

NB: no more used in the QE-app. We use phonopy instead.
"""

########################
################################ END DESCRIPTION
########################


def compute_pdos(
    fc,
):
    """Function to calculate square of a number.

    Args:
        fc: Force Constant object as obtained using the `generate_force_constant_instance` method.

    Returns:
        [data, parameters, type_of_data]: [pdos, parameters for plot,  type of data (bands, dos ...)].
    """
    # PDOS, we should allow changes in Ebins and mesh grid.
    # we should provide this as ProjectionData for the new BandsPlotWidget.

    # start Euphonic
    phonons, mode_grads = fc.calculate_qpoint_phonon_modes(
        util.mp_grid(
            [
                20,
                20,
                20,
            ]
        ),
        return_mode_gradients=True,
    )

    # adaptive broadening: https://euphonic.readthedocs.io/en/stable/dos.html#adaptive-broadening
    mode_widths = util.mode_gradients_to_widths(mode_grads, fc.crystal.cell_vectors)

    energy_bins = np.arange(0, 200, 0.2) * eu.ureg("meV")

    weightings = False

    if not weightings:
        pdos = phonons.calculate_pdos(energy_bins, mode_widths=mode_widths)

    else:
        pdos = phonons.calculate_pdos(
            energy_bins,
            mode_widths=mode_widths,
            weighting="coherent-plus-incoherent",
            cross_sections="BlueBook",
        )

    species_pdos = pdos.group_by("species")
    total_dos = pdos.sum()  # total dos

    # now we need to set the labels up properly
    for data in species_pdos.metadata["line_data"]:
        data["label"] = data["species"]
    # and then for the total
    total_dos.metadata["label"] = "Total"
    # end Euphonic

    dos = []
    # The total dos parsed
    tdos = {
        "label": "Total DOS",
        "x": total_dos.x_data.magnitude.tolist(),
        "y": total_dos.y_data.magnitude.tolist(),
        "borderColor": "#8A8A8A",  # dark gray
        "backgroundColor": "#8A8A8A",  # light gray
        "backgroundAlpha": "40%",
        "lineStyle": "solid",
    }
    dos.append(tdos)

    t = 0
    for atom in species_pdos.metadata["line_data"]:
        tdos = {
            "label": atom["label"],
            "x": total_dos.x_data.magnitude.tolist(),
            "y": species_pdos.y_data.magnitude[t].tolist(),
            "borderColor": cmap(atom["label"]),
            "backgroundColor": cmap(atom["label"]),
            "backgroundAlpha": "40%",
            "lineStyle": "solid",
        }
        t += 1
        dos.append(tdos)

    data_dict = {
        "fermi_energy": 0,  # I do not want it in my plot
        "dos": dos,
    }

    parameters = {}
    parameters["energy_range"] = {
        "ymin": np.min(total_dos.x_data.magnitude),
        "ymax": np.max(total_dos.x_data.magnitude),
    }

    return [data_dict, parameters]  # output_data.


def compute_bands(fc):
    """Function to calculate square of a number.

    Args:
        fc: Force Constant object as obtained using the `generate_force_constant_instance` method.

    Returns:
        [data, parameters, type_of_data]: [bands, parameters for plot,  type of data (bands, dos ...)].
    """

    # We need support for low-D systems.

    # start Euphonic
    cell = fc.crystal.to_spglib_cell()
    qpts = seekpath.get_explicit_k_path(cell)["explicit_kpoints_rel"]

    phonons = fc.calculate_qpoint_phonon_modes(qpts, asr="reciprocal")
    disp = phonons.get_dispersion()
    # end Euphonic

    bands = orm.BandsData()
    bands.set_kpoints(qpts)
    bands.set_bands(disp.y_data.magnitude.T, units="meV")

    curated_labels = []
    for label in disp.x_tick_labels:
        if label[1] == "$\\Gamma$":
            curated_labels.append((label[0], "Gamma"))
        else:
            curated_labels.append(label)

    bands.labels = curated_labels

    data = json.loads(bands._exportcontent("json", comments=False)[0])
    # The fermi energy from band calculation is not robust.
    """data["fermi_level"] = (
        fermi_energy or node.outputs.phonons.band_parameters["fermi_energy"]
    )"""
    # to be optimized: use the above results!!!
    data["fermi_level"] = 0
    data["Y_label"] = "Dispersion (meV)"

    # it does work now.
    parameters = {}

    bands = bands.get_bands()
    parameters["energy_range"] = {
        "ymin": np.min(bands) - 0.1,
        "ymax": np.max(bands) + 0.1,
    }

    return [jsanitize(data), parameters, "bands"]  # output_data.
