from __future__ import annotations

from aiidalab_qe.common.mvc import Model
import traitlets as tl
from aiida.common.extendeddicts import AttributeDict
from IPython.display import display
from euphonic import ForceConstants


class SingleCrystalFullModel(Model):
    node = tl.Instance(AttributeDict, allow_none=True)

    fc = tl.Instance(ForceConstants, allow_none=True)
    q_path = tl.Dict(allow_none=True)

    custom_kpath = tl.Unicode("")
    q_spacing = tl.Float(0.01)
    energy_broadening = tl.Float(0.05)
    energy_bins = tl.Int(200)
    temperature = tl.Float(0)
    weighting = tl.Unicode("coherent")

    E_units_button_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ("meV", "meV"),
            ("THz", "THz"),
        ],
    )
    E_units = tl.Unicode("meV")

    slider_intensity = tl.List(
        trait=tl.Float(),
        default_value=[1, 10],
    )

    parameters = tl.Dict()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.update_parameters()

        # Observe changes in dependent trailets
        self.observe(
            self.update_parameters,
            names=[
                "weighting",
                "E_units",
                "temperature",
                "q_spacing",
                "energy_broadening",
                "energy_bins",
            ],
        )

    def update_parameters(self):
        """Update the parameters dictionary dynamically."""
        self.parameters = {
            "weighting": self.weighting,
            "grid": None,
            "grid_spacing": 0.1,
            "energy_units": self.E_units,
            "temperature": self.temperature,
            "shape": "gauss",
            "length_unit": "angstrom",
            "q_spacing": self.q_spacing,
            "energy_broadening": self.energy_broadening,
            "q_broadening": None,
            "ebins": self.energy_bins,
            "e_min": 0,
            "e_max": None,
            "title": None,
            "ylabel": "THz",
            "xlabel": "",
            "save_json": False,
            "no_base_style": False,
            "style": False,
            "vmin": None,
            "vmax": None,
            "save_to": None,
            "asr": None,
            "dipole_parameter": 1.0,
            # "use_c": False,
            # "n_threads": None,
        }

    def _update_spectra(self):
        q_path = self.q_path
        if self.custom_kpath:
            q_path = self.q_path
            q_path["coordinates"], q_path["labels"] = self.curate_path_and_labels(
                self.custom_kpath
            )
            q_path["delta_q"] = self.q_spacing

        # spectra, parameters = produce_bands_weigthed_data(
        #     params=self.parameters,
        #     fc=self.fc,
        #     linear_path=q_path,
        #     plot=False,
        # )

        # if self.custom_path:
        #     self.x, self.y = np.meshgrid(
        #         spectra[0].x_data.magnitude, spectra[0].y_data.magnitude
        #     )
        #     (
        #         self.final_xspectra,
        #         self.final_zspectra,
        #         self.ticks_positions,
        #         self.ticks_labels,
        #     ) = generated_curated_data(spectra)
        # else:
        #     # Spectrum2D as output of the powder data
        #     self.x, self.y = np.meshgrid(
        #         spectra.x_data.magnitude, spectra.y_data.magnitude
        #     )

        #     # we don't need to curate the powder data,
        #     # we can directly use them:
        #     self.final_xspectra = spectra.x_data.magnitude
        #     self.final_zspectra = spectra.z_data.magnitude

    def curate_path_and_labels(self, path):
        # This is used to curate the path and labels of the spectra if custom kpath is provided.
        # I do not like this implementation (MB)
        coordinates = []
        labels = []
        linear_paths = path.split("|")
        for i in linear_paths:
            scoords = []
            s = i.split(
                " - "
            )  # not i.split("-"), otherwise also the minus of the negative numbers are used for the splitting.
            for k in s:
                labels.append(k.strip())
                # AAA missing support for fractions.
                l = tuple(map(float, [kk for kk in k.strip().split(" ")]))  # noqa: E741
                scoords.append(l)
            coordinates.append(scoords)
        return coordinates, labels

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
