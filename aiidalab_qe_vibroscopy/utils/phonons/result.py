"""Bands results view widgets

"""

from widget_bandsplot import BandsPlotWidget

from aiidalab_qe.common.panel import ResultPanel
from aiidalab_qe.plugins.pdos.result import cmap

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

        if "phonon_pdos" in node.outputs.vibronic:

            phonopy_calc = node.outputs.vibronic.phonon_pdos.creator

            kwargs = {}
            if "settings" in phonopy_calc.inputs:
                the_settings = phonopy_calc.inputs.settings.get_dict()
                for key in ["symmetrize_nac", "factor_nac", "subtract_residual_forces"]:
                    if key in the_settings:
                        kwargs.update({key: the_settings[key]})

            if "phonopy_data" in phonopy_calc.inputs:
                instance = phonopy_calc.inputs.phonopy_data.get_phonopy_instance(
                    **kwargs
                )

            elif "force_constants" in phonopy_calc.inputs:
                instance = phonopy_calc.inputs.force_constants.get_phonopy_instance(
                    **kwargs
                )

            symbols = instance.get_primitive().get_chemical_symbols()
            pdos = node.outputs.vibronic.phonon_pdos

            index_dict, dos_dict = {}, {
                "total_dos": np.zeros(np.shape(pdos.get_y()[0][1]))
            }
            for atom in set(symbols):
                # index lists
                index_dict[atom] = [
                    i for i in range(len(symbols)) if symbols[i] == atom
                ]
                # initialization of the pdos
                dos_dict[atom] = np.zeros(
                    np.shape(pdos.get_y()[index_dict[atom][0]][1])
                )

                for atom_contribution in index_dict[atom][:]:
                    dos_dict[atom] += pdos.get_y()[atom_contribution][1]
                    dos_dict["total_dos"] += pdos.get_y()[atom_contribution][1]

            dos = []
            # The total dos parsed
            tdos = {
                "label": "Total DOS",
                "x": pdos.get_x()[1].tolist(),
                "y": dos_dict.pop("total_dos").tolist(),
                "borderColor": "#8A8A8A",  # dark gray
                "backgroundColor": "#8A8A8A",  # light gray
                "backgroundAlpha": "40%",
                "lineStyle": "solid",
            }
            dos.append(tdos)

            t = 0
            for atom in dos_dict.keys():
                tdos = {
                    "label": atom,
                    "x": pdos.get_x()[1].tolist(),
                    "y": dos_dict[atom].tolist(),
                    "borderColor": cmap(atom),
                    "backgroundColor": cmap(atom),
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
                "ymin": np.min(dos[0]["x"]),
                "ymax": np.max(dos[0]["x"]),
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
