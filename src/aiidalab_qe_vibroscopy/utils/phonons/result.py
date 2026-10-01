"""Bands results view widgets"""

from aiidalab_qe.common.bands_pdos.utils import _cmap, _get_bands_labeling

import numpy as np
import json


def replace_symbols_with_uppercase(data):
    symbols_mapping = {
        "$\Gamma$": "\u0393",
        "$\\Gamma$": "\u0393",  # noqa: F601
        "$\\Delta$": "\u0394",
        "$\\Lambda$": "\u039b",
        "$\\Sigma$": "\u03a3",
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

    full_data = {
        "bands": None,
        "pdos": None,
        "thermo": None,
    }
    parameters = {}

    if "phonon_bands" in node.outputs:
        """
        copied and pasted from aiidalab_qe.common.bandsplotwidget.
        adapted for phonon outputs
        """

        data = json.loads(
            node.outputs.phonon_bands._exportcontent("json", comments=False)[0]
        )
        # The fermi energy from band calculation is not robust.
        data["fermi_energy"] = 0
        data["pathlabels"] = _get_bands_labeling(data)
        replace_symbols_with_uppercase(data["pathlabels"])
        data["Y_label"] = "Dispersion (THz)"

        bands = node.outputs.phonon_bands._get_bandplot_data(
            cartesian=True, prettify_format=None, join_symbol=None, get_segments=True
        )
        parameters["energy_range"] = {
            "ymin": np.min(bands["y"]) - 0.1,
            "ymax": np.max(bands["y"]) + 0.1,
        }
        data["band_type_idx"] = bands["band_type_idx"]
        data["x"] = bands["x"]
        data["y"] = bands["y"]
        full_data["bands"] = [data, parameters]

        if "phonon_pdos" in node.outputs:
            phonopy_calc = node.outputs.phonon_pdos.creator

            kwargs = {}
            if "settings" in phonopy_calc.inputs:
                the_settings = phonopy_calc.inputs.settings.get_dict()
                for key in ["symmetrize_nac", "factor_nac", "subtract_residual_forces"]:
                    if key in the_settings:
                        kwargs.update({key: the_settings[key]})

            symbols = node.inputs.structure.get_ase().get_chemical_symbols()
            pdos = node.outputs.phonon_pdos

            index_dict, dos_dict = (
                {},
                {"total_dos": np.zeros(np.shape(pdos.get_y()[0][1]))},
            )
            for i, atom in enumerate(set(symbols)):
                # index lists
                index_dict[atom] = [
                    j for j in range(len(symbols)) if symbols[j] == atom
                ]
                # initialization of the pdos
                dos_dict[atom] = np.zeros(np.shape(pdos.get_y()[i][1]))

                for atom_contribution in index_dict[atom][:]:
                    if len(pdos.get_y()) <= atom_contribution:
                        # I need this as for Al4, only one pdos tuple is provided...
                        # for Si2, actually, two are provided...
                        dos_dict[atom] += pdos.get_y()[i][1]
                        dos_dict["total_dos"] += pdos.get_y()[i][1]
                    else:
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
                    "borderColor": _cmap(atom),
                    "backgroundColor": _cmap(atom),
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

        if "phonon_thermo" in node.outputs:
            (
                what,
                T,
                units_k,
            ) = node.outputs.phonon_thermo.get_x()
            (
                F_name,
                F_data,
                units_F,
            ) = node.outputs.phonon_thermo.get_y()[0]
            (
                Entropy_name,
                Entropy_data,
                units_entropy,
            ) = node.outputs.phonon_thermo.get_y()[1]
            (
                Cv_name,
                Cv_data,
                units_Cv,
            ) = node.outputs.phonon_thermo.get_y()[2]

            full_data["thermo"] = (
                [T, F_data, units_F, Entropy_data, units_entropy, Cv_data, units_Cv],
                [],
                "thermal",
            )

        return full_data
    else:
        return None


def atomic_pdos_from_workchain(node):
    """Read per-atom PDOS in displayed unit-cell order without recalculation."""
    from aiidalab_qe_vibroscopy.utils.mode_projection import unitcell_to_primitive

    pdos = node.outputs.phonon_pdos
    creator = pdos.creator
    if creator is None:
        raise ValueError("The source of the saved atom-resolved DOS is unavailable.")
    parameters = {
        key.lower(): value
        for key, value in creator.inputs.parameters.get_dict().items()
    }
    if str(parameters.get("pdos", "")).lower() != "auto":
        raise ValueError(
            "This calculation did not save DOS in automatic per-atom order."
        )
    if "primitive_axes" in parameters:
        raise ValueError(
            "Atom mapping for a manually overridden primitive cell is unavailable."
        )
    source = next(
        (
            getattr(creator.inputs, key)
            for key in ("phonopy_data", "force_constants")
            if key in creator.inputs
        ),
        None,
    )
    if source is None or not hasattr(source, "get_phonopy_instance"):
        raise ValueError(
            "The stored phonon data needed for atom mapping are unavailable."
        )
    phonopy = source.get_phonopy_instance()
    structure = node.inputs.structure.get_ase()
    mapping = unitcell_to_primitive(phonopy, structure)
    values = np.array([row[1] for row in pdos.get_y()])
    number_of_primitive_atoms = len(phonopy.primitive)
    if (
        parameters.get("xyz_projection")
        and len(values) == 3 * number_of_primitive_atoms
    ):
        values = values.reshape(number_of_primitive_atoms, 3, -1).sum(axis=1)
    if len(values) != number_of_primitive_atoms:
        raise ValueError("The saved DOS curves cannot be assigned to individual atoms.")
    # Report DOS per displayed unit cell, as in the existing element-grouped view.
    return values[mapping]
