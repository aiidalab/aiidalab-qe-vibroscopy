"""Bands results view widgets

"""


from aiidalab_qe.common.panel import ResultPanel
from aiidalab_qe.common.bandpdoswidget import cmap, get_bands_labeling

import numpy as np
import json


def replace_symbols_with_uppercase(data):
    symbols_mapping = {
        "$\Gamma$": "\u0393",
        "$\\Gamma$": "\u0393",
        "$\\Delta$": "\u0394",
        "$\\Lambda$": "\u039B",
        "$\\Sigma$": "\u03A3",
        "$\\Epsilon$": "\u0395",
    }

    for sublist in data:
        for i, element in enumerate(sublist):
            if element in symbols_mapping:
                sublist[i] = symbols_mapping[element]


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

        """
        copied and pasted from aiidalab_qe.common.bandsplotwidget.
        adapted for phonon outputs
        """

        data = json.loads(
            node.outputs.vibronic.phonon_bands._exportcontent("json", comments=False)[0]
        )
        # The fermi energy from band calculation is not robust.
        data["fermi_energy"] = 0
        data["pathlabels"] = get_bands_labeling(data)
        replace_symbols_with_uppercase(data["pathlabels"])
        data["Y_label"] = "Dispersion (THz)"

        bands = node.outputs.vibronic.phonon_bands._get_bandplot_data(cartesian=True, prettify_format=None, join_symbol=None, get_segments=True)
        parameters["energy_range"] = {
            "ymin": np.min(bands["y"]) - 0.1,
            "ymax": np.max(bands["y"]) + 0.1,
        }
        data["band_type_idx"] = bands["band_type_idx"]
        data["x"] = bands["x"]
        data["y"] = bands["y"]
        full_data["bands"] = [data, parameters]

        if "phonon_pdos" in node.outputs.vibronic:

            phonopy_calc = node.outputs.vibronic.phonon_pdos.creator

            kwargs = {}
            if "settings" in phonopy_calc.inputs:
                the_settings = phonopy_calc.inputs.settings.get_dict()
                for key in ["symmetrize_nac", "factor_nac", "subtract_residual_forces"]:
                    if key in the_settings:
                        kwargs.update({key: the_settings[key]})

            symbols = node.inputs.structure.get_ase().get_chemical_symbols()
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
