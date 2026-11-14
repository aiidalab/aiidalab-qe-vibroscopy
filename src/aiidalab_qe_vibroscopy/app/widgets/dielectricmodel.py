from aiidalab_qe.common.mvc import Model
import traitlets as tl
from aiida.common.extendeddicts import AttributeDict
from aiidalab_qe_vibroscopy.utils.dielectric.result import NumpyEncoder
import numpy as np
import base64
import json
from IPython.display import display
from aiidalab_qe_vibroscopy.utils.dielectric.result import export_dielectric_data


class DielectricModel(Model):
    vibro = tl.Instance(AttributeDict, allow_none=True)

    site_selector_options = tl.List(
        trait=tl.Tuple((tl.Unicode(), tl.Int())),
    )
    dielectric_data = {}

    dielectric_tensor_table = tl.Unicode("")
    born_charges_table = tl.Unicode("")
    raman_tensors_table = tl.Unicode("")
    site = tl.Int()

    def fetch_data(self):
        """Fetch the dielectric data from the VibroWorkChain"""
        self.dielectric_data = export_dielectric_data(self.vibro)

    def set_initial_values(self):
        """Set the initial values for the model."""

        self.dielectric_tensor_table = self._create_dielectric_tensor_table()
        self.born_charges_table = self._create_born_charges_table(0)
        self.raman_tensors_table = self._create_raman_tensors_table(0)
        self.site_selector_options = self._get_site_selector_options()

    def _get_site_selector_options(self):
        """Get the site selector options."""
        if not self.dielectric_data:
            return []

        unit_cell_sites = self.dielectric_data["unit_cell"]
        decimal_places = 5
        # Create the options with rounded positions
        site_selector_options = [
            (
                f"{site.kind_name} @ ({', '.join(f'{coord:.{decimal_places}f}' for coord in site.position)})",
                index,
            )
            for index, site in enumerate(unit_cell_sites)
        ]
        return site_selector_options

    def _create_dielectric_tensor_table(self):
        """Create the HTML table for the dielectric tensor."""
        if not self.dielectric_data:
            return ""

        dielectric_tensor = self.dielectric_data["dielectric_tensor"]
        table_data = self._generate_table(dielectric_tensor)
        return table_data

    def _create_born_charges_table(self, site_index):
        """Create the HTML table for the Born charges."""
        if not self.dielectric_data:
            return ""

        born_charges = self.dielectric_data["born_charges"]
        round_data = born_charges[site_index].round(6)
        table_data = self._generate_table(round_data)
        return table_data

    def _create_raman_tensors_table(self, site_index):
        """Create the HTML table for the Raman tensors."""
        if not self.dielectric_data:
            return ""

        raman_tensors = self.dielectric_data["raman_tensors"]
        round_data = raman_tensors[site_index].round(6)
        table_data = self._generate_table(round_data, cell_width="200px")
        return table_data

    def download_data(self, _=None):
        """Function to download the data."""
        if self.dielectric_data:
            data_to_print = {
                key: value
                for key, value in self.dielectric_data.items()
                if key != "unit_cell"
            }
            file_name = "dielectric_data.json"
            json_str = json.dumps(data_to_print, cls=NumpyEncoder)
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

    def on_site_selection_change(self, site):
        self.site = site
        self.born_charges_table = self._create_born_charges_table(site)
        self.raman_tensors_table = self._create_raman_tensors_table(site)

    def _generate_table(self, data, cell_width="50px"):
        rows = []
        for row in data:
            cells = []
            for value in row:
                # Check if value is a numpy array
                if isinstance(value, np.ndarray):
                    # Format the numpy array as a string, e.g., "[0, 0, 1]"
                    value_str = np.array2string(
                        value, separator=", ", formatter={"all": lambda x: f"{x:.6g}"}
                    )
                    cell = f"<td>{value_str}</td>"
                elif isinstance(value, str) and value == "special":
                    # Handle the "special" keyword
                    cell = f'<td class="blue-cell">{value}</td>'
                else:
                    # Handle other types (numbers, strings, etc.)
                    cell = f"<td>{value}</td>"
                cells.append(cell)
            rows.append(f"<tr>{''.join(cells)}</tr>")

        # Define the HTML with styles, using the dynamic cell width
        table_html = f"""
        <style>
            table {{
                border-collapse: collapse;
                width: auto; /* Adjust to content */
            }}
            th, td {{
                border: 1px solid black;
                text-align: center;
                padding: 4px;
                height: 12px;
                width: {cell_width}; /* Set custom cell width */
            }}
            .blue-cell {{
                background-color: gray;
                color: white;
            }}
        </style>
        <table>
            {"".join(rows)}
        </table>
        """
        return table_html
