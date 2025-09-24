# -*- coding: utf-8 -*-

import json

import numpy as np


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()  # Convert ndarray to list
        return super().default(obj)


def create_html_table(matrix):
    """
    Create an HTML table representation of a 3x3 matrix.

    :param matrix: List of lists representing a 3x3 matrix
    :return: HTML table string
    """
    html = '<table border="1" style="border-collapse: collapse;">'
    for row in matrix:
        html += "<tr>"
        for cell in row:
            html += f'<td style="padding: 5px; text-align: center;">{cell}</td>'
        html += "</tr>"
    html += "</table>"
    return html


def get_priority_tensor(filtered_node):
    """
    Retrieve the tensor from filtered_nodes based on the predefined priority of keys.

    :param filtered_node: Node with outputs containing tensors or VibrationalData
    :return: Corresponding to the highest priority key found, or None if not found
    """
    # Define the priority order of keys within the function
    priority_keys = [
        "numerical_accuracy_4",
        "numerical_accuracy_2_step_2",
        "numerical_accuracy_2_step_1",
        "numerical_accuracy_2",
    ]

    # Get the keys from the tensor outputs
    tensor_keys = filtered_node.keys()

    for key in priority_keys:
        if key in tensor_keys:
            return filtered_node[key]

    # If no matching key is found, return None or handle the case as needed
    return None


def export_dielectric_data(node):
    if not any(key in node for key in ["iraman", "dielectric", "harmonic"]):
        return None

    else:
        if "iraman" in node:
            vibrational_data = node.iraman.vibrational_data

        elif "harmonic" in node:
            vibrational_data = node.harmonic.vibrational_data

        elif "dielectric" in node:
            tensor_data = node.dielectric
            output_data = get_priority_tensor(getattr(tensor_data, "tensors", tensor_data))
            dielectric_tensor = output_data.get_array("dielectric").round(
                6
            )  # Dielectric Constant
            born_charges = output_data.get_array(
                "born_charges"
            )  # List of Born effective charges per Atom
            vol = output_data.get_unitcell().get_cell_volume()  # Volume of the cell
            raman_tensors = output_data.get_array("raman_tensors")  # Raman tensors
            nlo_susceptibility = output_data.get_array(
                "nlo_susceptibility"
            )  # non-linear optical susceptibility tensor (pm/V)
            unit_cell = output_data.get_unitcell().sites
            return {
                "dielectric_tensor": dielectric_tensor,
                "born_charges": born_charges,
                "volume": vol,
                "raman_tensors": raman_tensors,
                "nlo_susceptibility": nlo_susceptibility,
                "unit_cell": unit_cell,
            }

        output_data = get_priority_tensor(vibrational_data)
        dielectric_tensor = output_data.dielectric.round(6)  # Dielectric Constant
        born_charges = output_data.get_array(
            "born_charges"
        )  # List of Born effective charges per Atom
        vol = output_data.get_unitcell().get_cell_volume()  # Volume of the cell
        raman_tensors = output_data.get_array("raman_tensors")  # Raman tensors
        nlo_susceptibility = (
            output_data.nlo_susceptibility
        )  # non-linear optical susceptibility tensor (pm/V)
        unit_cell = output_data.get_unitcell().sites

        return {
            "dielectric_tensor": dielectric_tensor,
            "born_charges": born_charges,
            "volume": vol,
            "raman_tensors": raman_tensors,
            "nlo_susceptibility": nlo_susceptibility,
            "unit_cell": unit_cell,
        }
