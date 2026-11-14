from aiidalab_qe.common.mvc import Model
import traitlets as tl

from aiida.orm.nodes.process.workflow.workchain import WorkChainNode
from aiidalab_qe_vibroscopy.utils.phonons.result import export_phononworkchain_data
from IPython.display import display
import numpy as np


class PhononModel(Model):
    vibro = tl.Instance(WorkChainNode, allow_none=True)

    pdos_data = {}
    bands_data = {}
    thermo_data = {}

    def fetch_data(self):
        """Fetch the phonon data from the VibroWorkChain"""
        phonon_data = export_phononworkchain_data(self.vibro)
        self.pdos_data = phonon_data["pdos"][0]
        self.bands_data = phonon_data["bands"][0]
        self.thermo_data = phonon_data["thermo"][0]

    def update_thermo_plot(self, fig):
        """Update the thermal properties plot."""
        self.temperature = self.thermo_data[0]
        self.free_E = self.thermo_data[1]
        F_units = self.thermo_data[2]
        self.entropy = self.thermo_data[3]
        E_units = self.thermo_data[4]
        self.Cv = self.thermo_data[5]
        Cv_units = self.thermo_data[6]
        fig.update_layout(
            xaxis=dict(
                title="Temperature (K)",
                linecolor="black",
                linewidth=2,
                showline=True,
            ),
            yaxis=dict(linecolor="black", linewidth=2, showline=True),
            plot_bgcolor="white",
        )
        fig.add_scatter(
            x=self.temperature, y=self.free_E, name=f"Helmoltz Free Energy ({F_units})"
        )
        fig.add_scatter(x=self.temperature, y=self.entropy, name=f"Entropy ({E_units})")
        fig.add_scatter(
            x=self.temperature, y=self.Cv, name=f"Specific Heat-V=const ({Cv_units})"
        )

    def download_thermo_data(self, _=None):
        """Function to download the phonon data."""
        import json
        import base64

        file_name = "phonon_thermo_data.json"
        data_export = {}
        for key, value in zip(
            [
                "Temperature (K)",
                "Helmoltz Free Energy (kJ/mol)",
                "Entropy (J/K/mol)",
                "Specific Heat-V=const (J/K/mol)",
            ],
            [self.temperature, self.free_E, self.entropy, self.Cv],
        ):
            if isinstance(value, np.ndarray):
                data_export[key] = value.tolist()
            else:
                data_export[key] = value

        json_str = json.dumps(data_export)
        b64_str = base64.b64encode(json_str.encode()).decode()
        self._download(payload=b64_str, filename=file_name)

    def download_bandspdos_data(self, _=None):
        """Function to download the phonon data."""
        import json
        from monty.json import jsanitize
        import base64

        file_name_bands = "phonon_bands_data.json"
        file_name_pdos = "phonon_dos_data.json"
        if self.bands_data:
            bands_data_export = {}
            for key, value in self.bands_data.items():
                if isinstance(value, np.ndarray):
                    bands_data_export[key] = value.tolist()
                else:
                    bands_data_export[key] = value

            json_str = json.dumps(jsanitize(bands_data_export))
            b64_str = base64.b64encode(json_str.encode()).decode()
            self._download(payload=b64_str, filename=file_name_bands)
        if self.pdos_data:
            json_str = json.dumps(jsanitize(self.pdos_data))
            b64_str = base64.b64encode(json_str.encode()).decode()
            self._download(payload=b64_str, filename=file_name_pdos)

    @staticmethod
    def _download(payload, filename):
        from IPython.display import Javascript

        javas = Javascript(
            """
            var link = document.createElement('a');
            link.href = 'data:text/json;charset=utf-8;base64,{payload}'
            link.download = "{filename}"
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            """.format(payload=payload, filename=filename)
        )
        display(javas)
