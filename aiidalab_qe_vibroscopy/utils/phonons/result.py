"""Bands results view widgets

"""

from widget_bandsplot import BandsPlotWidget

from aiidalab_qe.common.panel import ResultPanel

import numpy as np


def export_phononworkchain_data(node, fermi_energy=None):

    """
    We have multiple choices: BANDS, DOS, THERMODYNAMIC.
    """

    import json

    from monty.json import jsanitize

    full_data = {
        "bands": None,
        "pdos": None,
        "thermo": None,
    }
    parameters = {}

    if not "vibronic" in node.outputs:
        return None

    if "phonon_bands" in node.outputs.vibronic:

        data = json.loads(
            node.outputs.vibronic.phonon_bands._exportcontent("json", comments=False)[0]
        )
        # The fermi energy from band calculation is not robust.
        """data["fermi_level"] = (
            fermi_energy or node.outputs.phonons.band_parameters["fermi_energy"]
        )"""
        # to be optimized: use the above results!!!
        bands = node.outputs.vibronic.phonon_bands.get_bands()
        data["fermi_level"] = 0
        data["Y_label"] = "Dispersion (THz)"

        # it does work now.
        parameters["energy_range"] = {
            "ymin": np.min(bands) - 0.1,
            "ymax": np.max(bands) + 0.1,
        }

        # TODO: THERMOD, FORCES; minors: bands-labels, done: no-fermi-in-dos.

        full_data["bands"] = [jsanitize(data), parameters]

        if "total_phonon_dos" in node.outputs.vibronic:
            (
                what,
                energy_dos,
                units_omega,
            ) = node.outputs.vibronic.phonon_dos.get_x()
            (
                dos_name,
                dos_data,
                units_dos,
            ) = node.outputs.vibronic.phonon_dos.get_y()[0]
            dos = []
            # The total dos parsed
            tdos = {
                "label": "Total DOS",
                "x": energy_dos.tolist(),
                "y": dos_data.tolist(),
                "borderColor": "#8A8A8A",  # dark gray
                "backgroundColor": "#8A8A8A",  # light gray
                "backgroundAlpha": "40%",
                "lineStyle": "solid",
            }
            dos.append(tdos)

            parameters["energy_range"] = {
                "ymin": np.min(energy_dos) - 0.1,
                "ymax": np.max(energy_dos) + 0.1,
            }

            data_dict = {
                "fermi_energy": 0,  # I do not want it in my plot
                "dos": dos,
            }

            full_data["pdos"] = [json.loads(json.dumps(data_dict)), parameters, "dos"]

        if "phonon_thermo" in node.outputs.vibronic:
            (
                what,
                T,
                units_k,
            ) = node.outputs.vibronic.phonon_thermo.get_x()
            (
                F_name,
                F_data,
                units_F,
            ) = node.outputs.vibronic.phonon_thermo.get_y()[0]
            (
                Entropy_name,
                Entropy_data,
                units_entropy,
            ) = node.outputs.vibronic.phonon_thermo.get_y()[1]
            (
                Cv_name,
                Cv_data,
                units_Cv,
            ) = node.outputs.vibronic.phonon_thermo.get_y()[2]

            full_data["thermo"] = (
                [T, F_data, units_F, Entropy_data, units_entropy, Cv_data, units_Cv],
                [],
                "thermal",
            )

        return full_data
    else:
        return None
