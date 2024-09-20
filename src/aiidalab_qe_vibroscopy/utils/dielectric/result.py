# -*- coding: utf-8 -*-

import ipywidgets as ipw
from IPython.display import HTML, display
import base64
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
    if "vibronic" not in node.outputs:
        return None

    if not any(
        key in node.outputs.vibronic for key in ["iraman", "dielectric", "harmonic"]
    ):
        return None

    else:
        if "iraman" in node.outputs.vibronic:
            vibrational_data = node.outputs.vibronic.iraman.vibrational_data

        elif "harmonic" in node.outputs.vibronic:
            vibrational_data = node.outputs.vibronic.harmonic.vibrational_data

        elif "dielectric" in node.outputs.vibronic:
            tensor_data = node.outputs.vibronic.dielectric
            output_data = get_priority_tensor(tensor_data)
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


class DielectricResults(ipw.VBox):
    def __init__(self, dielectric_data):
        # Helper
        self.dielectric_results_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            The DielectricWorkchain computes different properties: <br>
                <em style="display: inline-block; margin-left: 20px;">-High Freq. Dielectric Tensor </em> <br>
                <em style="display: inline-block; margin-left: 20px;">-Born Charges </em> <br>
                <em style="display: inline-block; margin-left: 20px;">-Raman Tensors </em> <br>
                <em style="display: inline-block; margin-left: 20px;">-The non-linear optical susceptibility tensor </em> <br>
                All information can be downloaded as a JSON file. <br>

            </div>"""
        )

        self.unit_cell_sites = dielectric_data.pop("unit_cell")
        self.dielectric_data = dielectric_data
        self.dielectric_tensor = dielectric_data["dielectric_tensor"]
        self.born_charges = dielectric_data["born_charges"]
        self.volume = dielectric_data["volume"]
        self.raman_tensors = dielectric_data["raman_tensors"]
        self.nlo_susceptibility = dielectric_data["nlo_susceptibility"]

        # HTML table with the dielectric tensor
        self.dielectric_tensor_table = ipw.Output()

        # HTML table with the Born charges @ site
        self.born_charges_table = ipw.Output()

        # HTML table with the Raman tensors @ site
        self.raman_tensors_table = ipw.Output()

        decimal_places = 6
        # Create the options with rounded positions
        site_selector_options = [
            (
                f"{site.kind_name} @ ({', '.join(f'{coord:.{decimal_places}f}' for coord in site.position)})",
                index,
            )
            for index, site in enumerate(self.unit_cell_sites)
        ]

        self.site_selector = ipw.Dropdown(
            options=site_selector_options,
            value=site_selector_options[0][1],
            layout=ipw.Layout(width="450px"),
            description="Select atom site:",
            style={"description_width": "initial"},
        )
        # Download button
        self.download_button = ipw.Button(
            description="Download Data", icon="download", button_style="primary"
        )
        self.download_button.on_click(self.download_data)

        # Initialize the HTML table
        self._create_dielectric_tensor_table()
        # Initialize Born Charges Table
        self._create_born_charges_table(self.site_selector.value)
        # Initialize Raman Tensors Table
        self._create_raman_tensors_table(self.site_selector.value)

        self.site_selector.observe(self._on_site_selection_change, names="value")
        super().__init__(
            children=(
                self.dielectric_results_help,
                ipw.HTML("<h3>Dielectric tensor</h3>"),
                self.dielectric_tensor_table,
                self.site_selector,
                ipw.HBox(
                    [
                        ipw.VBox(
                            [
                                ipw.HTML("<h3>Born effective charges</h3>"),
                                self.born_charges_table,
                            ]
                        ),
                        ipw.VBox(
                            [
                                ipw.HTML("<h3>Raman Tensor </h3>"),
                                self.raman_tensors_table,
                            ]
                        ),
                    ]
                ),
                self.download_button,
            )
        )

    def download_data(self, _=None):
        """Function to download the data."""
        file_name = "dielectric_data.json"

        json_str = json.dumps(self.dielectric_data, cls=NumpyEncoder)
        b64_str = base64.b64encode(json_str.encode()).decode()
        self._download(payload=b64_str, filename=file_name)

    @staticmethod
    def _download(payload, filename):
        """Download payload as a file named as filename."""
        from IPython.display import Javascript

        javas = Javascript(
            f"""
            var link = document.createElement('a');
            link.href = 'data:text/json;charset=utf-8;base64,{payload}'
            link.download = "{filename}"
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            """
        )
        display(javas)

    def _create_html_table(self, matrix):
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

    def _create_dielectric_tensor_table(self):
        table_data = self._create_html_table(self.dielectric_tensor)
        self.dielectric_tensor_table.layout = {
            "overflow": "auto",
            "height": "100px",
            "width": "300px",
        }
        with self.dielectric_tensor_table:
            display(HTML(table_data))

    def _create_born_charges_table(self, site_index):
        round_data = self.born_charges[site_index].round(6)
        table_data = self._create_html_table(round_data)
        self.born_charges_table.layout = {
            "overflow": "auto",
            "height": "150px",
            "width": "300px",
        }
        with self.born_charges_table:
            display(HTML(table_data))

    def _create_raman_tensors_table(self, site_index):
        round_data = self.raman_tensors[site_index].round(6)
        table_data = self._create_html_table(round_data)
        self.raman_tensors_table.layout = {
            "overflow": "auto",
            "height": "150px",
            "width": "500px",
        }
        with self.raman_tensors_table:
            display(HTML(table_data))

    def _on_site_selection_change(self, change):
        self.born_charges_table.clear_output()
        self.raman_tensors_table.clear_output()
        self._create_born_charges_table(change["new"])
        self._create_raman_tensors_table(change["new"])
