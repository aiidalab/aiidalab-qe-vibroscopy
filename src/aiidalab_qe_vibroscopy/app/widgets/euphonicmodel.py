import numpy as np
import traitlets as tl
import copy
from IPython.display import display

from aiidalab_qe_vibroscopy.utils.euphonic.data.structure_factors import (
    AttrDict,
    produce_bands_weigthed_data,
    produce_powder_data,
    generated_curated_data,
    produce_Q_section_spectrum,
    produce_Q_section_modes,
)

from aiidalab_qe_vibroscopy.utils.euphonic.data.phonopy_interface import (
    generate_force_constant_from_phonopy,
)

from aiidalab_qe_vibroscopy.utils.euphonic.data.export_vibronic_to_euphonic import (
    export_euphonic_data,
)

from aiidalab_qe_vibroscopy.utils.euphonic.data.parameters import (
    parameters_single_crystal,
    parameters_powder,
)


from aiidalab_qe.common.mvc import Model


class EuphonicResultsModel(Model):
    """Model for the neutron scattering results panel."""

    # Here we mode all the model and data-controller, i.e. all the data and their
    # manipulation to produce new spectra.
    # plot-controller, i.e. scales, colors and so on, should be attached to the widget, I think.
    # the above last point should be discussed more.

    # For now, we do only the first of the following:
    # 1. single crystal data: sc
    # 2. powder average: pa
    # 3. Q planes: qp

    # Settings for single crystal and powder average
    q_spacing = tl.Float(0.1)  # q-spacing for the linear path
    energy_broadening = tl.Float(0.5)  # energy broadening
    energy_bins = tl.Int(200)  # energy bins
    temperature = tl.Float(0)  # temperature
    weighting = tl.Unicode("coherent")  # weighting
    energy_units = tl.Unicode("meV")  # energy units
    intensity_filter = tl.List(
        trait=tl.Float(), default_value=[0, 100]
    )  # intensity filter

    THz_to_meV = 4.13566553853599  # conversion factor.
    THz_to_cm1 = 33.3564095198155  # conversion factor.

    def __init__(
        self,
        node=None,
        spectrum_type: str = "single_crystal",
        detached_app: bool = False,
        **kwargs,
    ):
        super().__init__(**kwargs)

        self.spectra = {}
        self.path = []
        self.q_path = None
        self.spectrum_type = spectrum_type
        self.xlabel = None
        self.ylabel = self.energy_units
        self.detached_app = detached_app
        if node:  # qe app mode.
            self.vibro = node

        if self.spectrum_type == "single_crystal":
            self._inject_single_crystal_settings()
        elif self.spectrum_type == "powder":
            self._inject_powder_settings()
        elif self.spectrum_type == "q_planes":
            self._inject_qsection_settings()

    def set_model_state(self, parameters: dict):
        for k, v in parameters.items():
            setattr(self, k, v)

    def _get_default(self, trait):
        if trait in ["h_vec", "k_vec"]:
            return [1, 1, 1, 100, 1]
        elif trait == "Q0_vec":
            return [0.0, 0.0, 0.0]
        elif trait == "intensity_filter":
            return [0, 100]
        return self.traits()[trait].default_value

    def get_model_state(self):
        return {trait: getattr(self, trait) for trait in self.traits()}

    def reset(
        self,
    ):
        with self.hold_trait_notifications():
            for trait in self.traits():
                if trait not in ["intensity_filter", "energy_units"]:
                    setattr(self, trait, self._get_default(trait))

    def fetch_data(self):
        """Fetch the data from the database or from the uploaded files."""
        # 1. from aiida, so we have the node
        if hasattr(self, "fc"):
            # we already have the data (this happens if I clone the model with already the data inside)
            return
        if self.vibro:
            ins_data = export_euphonic_data(self.vibro)
            self.fc = ins_data["fc"]
            self.q_path = ins_data["q_path"]
        # 2. from uploaded files...
        else:
            # here we just use upload_widget as MVC all together, for simplicity.
            # moreover, this part is not used in the current QE app.
            self.fc = self.upload_widget._read_phonopy_files(
                fname=self.fname,
                phonopy_yaml_content=self.phonopy_yaml_content,
                fc_hdf5_content=self.fc_hdf5_content,
            )

    def _inject_single_crystal_settings(
        self,
    ):
        # Case in which we want to inject the model into the single crystal widget
        # we define specific parameters dictionary and callback function for the single crystal case
        self.parameters = copy.deepcopy(parameters_single_crystal)
        self._callback_spectra_generation = produce_bands_weigthed_data

        # Dynamically add a trait for single crystal settings
        self.add_traits(custom_kpath=tl.Unicode(""))

    def _inject_powder_settings(
        self,
    ):
        self.parameters = copy.deepcopy(parameters_powder)
        self._callback_spectra_generation = produce_powder_data

        # Dynamically add a trait for powder settings
        self.add_traits(q_min=tl.Float(0.0))
        self.add_traits(q_max=tl.Float(1))
        self.add_traits(npts=tl.Int(500))

    def _inject_qsection_settings(
        self,
    ):
        # self._callback_spectra_generation = produce_Q_section_modes
        # Dynamically add a trait for q section settings
        self.add_traits(center_e=tl.Float(0.0))
        self.add_traits(Q0_vec=tl.List(trait=tl.Float(), default_value=[0.0, 0.0, 0.0]))
        self.add_traits(
            h_vec=tl.List(trait=tl.Float(), default_value=[1, 1, 1, 100, 1])
        )
        self.add_traits(
            k_vec=tl.List(trait=tl.Float(), default_value=[1, 1, 1, 100, 1])
        )

    def get_spectra(
        self,
    ):
        # This is used to update the spectra when the parameters are changed

        if self.spectrum_type == "q_planes":
            self._get_qsection_spectra()
        else:
            self.parameters.update(self.get_model_state())
            # custom linear path
            custom_kpath = self.custom_kpath if hasattr(self, "custom_kpath") else ""
            if len(custom_kpath) > 1:
                coordinates, labels = self._curate_path_and_labels()
                qpath = {
                    "coordinates": coordinates,
                    "labels": labels,  # ["$\Gamma$","X","X","(1,1,1)"],
                    "delta_q": self.parameters["q_spacing"],
                }
            else:
                qpath = copy.deepcopy(self.q_path)
                if qpath:
                    qpath["delta_q"] = self.parameters["q_spacing"]

            spectra, parameters = self._callback_spectra_generation(
                params=AttrDict(self.parameters),
                fc=self.fc,
                linear_path=qpath,
                plot=False,
            )

        # curated spectra (labels and so on...)
        if self.spectrum_type == "single_crystal":  # single crystal case
            self.x, self.y = np.meshgrid(
                spectra.x_data.magnitude, spectra.y_data.magnitude
            )
            (
                final_xspectra,
                final_zspectra,
                ticks_positions,
                ticks_labels,
            ) = generated_curated_data(spectra)

            self.ticks_positions = ticks_positions
            self.ticks_labels = ticks_labels

            self.z = final_zspectra.T
            self.y = self.y[:, 0]
            self.x = list(
                range(self.ticks_positions[-1] + 1)
            )  # we have, instead, the ticks positions and labels

            # we need to cut out some of the x and y data, as they are not used in the plot.
            self.y = self.y[: np.shape(self.z)[0]]
            self.x = self.x[: np.shape(self.z)[1]]

            self.ylabel = self.energy_units

        elif self.spectrum_type == "powder":  # powder case
            # Spectrum2D as output of the powder data
            self.x, self.y = np.meshgrid(
                spectra.x_data.magnitude, spectra.y_data.magnitude
            )

            # we don't need to curate the powder data, at variance with the single crystal case.
            # We can directly use them:
            self.xlabel = "|q| (1/A)"
            self.x = spectra.x_data.magnitude
            self.y = self.y[:, 0]
            self.z = spectra.z_data.magnitude.T

            # we need to cut out some of the x and y data, as they are not used in the plot.
            self.y = self.y[: np.shape(self.z)[0]]
            self.x = self.x[: np.shape(self.z)[1]]

        elif self.spectrum_type == "q_planes":
            pass
        else:
            raise ValueError("Spectrum type not recognized:", self.spectrum_type)

        self.y = self.y * self.energy_conversion_factor(self.energy_units, "meV")

    def _get_qsection_spectra(
        self,
    ):
        # This is used to update the spectra in the case we plot the Q planes (the third tab).
        self.parameters_qplanes = AttrDict(
            {
                "h": np.array([i for i in self.h_vec[:-2]]),
                "k": np.array([i for i in self.k_vec[:-2]]),
                "n_h": int(self.h_vec[-2]),
                "n_k": int(self.k_vec[-2]),
                "h_extension": self.h_vec[-1],
                "k_extension": self.k_vec[-1],
                "Q0": np.array([i for i in self.Q0_vec[:]]),
                "ecenter": self.center_e,
                "deltaE": self.energy_broadening,
                "bins": self.energy_bins,
                "spectrum_type": self.weighting,
                "temperature": self.temperature,
            }
        )

        modes, q_array, h_array, k_array, labels, dw = produce_Q_section_modes(
            self.fc,
            h=self.parameters_qplanes.h,
            k=self.parameters_qplanes.k,
            Q0=self.parameters_qplanes.Q0,
            n_h=self.parameters_qplanes.n_h,
            n_k=self.parameters_qplanes.n_k,
            h_extension=self.parameters_qplanes.h_extension,
            k_extension=self.parameters_qplanes.k_extension,
            temperature=self.parameters_qplanes.temperature,
        )

        self.z, q_array, self.x, self.y, self.labels = produce_Q_section_spectrum(
            modes,
            q_array,
            h_array,
            k_array,
            ecenter=self.parameters_qplanes.ecenter,
            deltaE=self.parameters_qplanes.deltaE,
            bins=self.parameters_qplanes.bins,
            spectrum_type=self.parameters_qplanes.spectrum_type,
            dw=dw,
            labels=labels,
        )
        self.xlabel = self.labels["h"]
        self.ylabel = self.labels["k"]

    def energy_conversion_factor(self, new, old):
        # TODO: check this is correct.
        if new == old:
            return 1
        if new == "meV":
            if old == "THz":
                return self.THz_to_meV
            elif old == "1/cm":
                return 1 / self.THz_to_cm1 * self.THz_to_meV
        elif new == "THz":
            if old == "meV":
                return 1 / self.THz_to_meV
            elif old == "1/cm":
                return 1 / self.THz_to_cm1
        elif new == "1/cm":
            if old == "meV":
                return 1 / self.THz_to_meV * self.THz_to_cm1
            elif old == "THz":
                return self.THz_to_cm1

    def _update_energy_units(self, new, old):
        # This is used to update the energy units in the plot.
        self.y = self.y * self.energy_conversion_factor(new, old)
        self.ylabel = self.energy_units

    def _curate_path_and_labels(
        self,
    ):
        # This is used to curate the path and labels of the spectra if custom kpath is provided.
        # I do not like this implementation (MB)
        coordinates = []
        labels = []
        path = self.custom_kpath
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

    def produce_phonopy_files(self):
        # This is used to produce the phonopy files from
        # PhonopyCalculation data. The files are phonopy.yaml and force_constants.hdf5
        phonopy_yaml, fc_hdf5 = generate_force_constant_from_phonopy(
            self.vibro.phonon_bands.creator, mode="download"
        )
        return phonopy_yaml, fc_hdf5

    def prepare_data_for_download(self):
        import pandas as pd
        import base64
        from aiidalab_qe_vibroscopy.utils.euphonic.plotting.generator import (
            generate_from_template,
        )

        random_number = np.random.randint(0, 100)

        # Plotted_data
        if self.spectrum_type == "q_planes":
            # we store x,y,z as values, not as values, indexes and columns
            z = self.z.reshape(int(self.k_vec[-2]) + 1, int(self.h_vec[-2]) + 1)
            df = pd.DataFrame(
                z,
                index=self.y[: int(self.k_vec[-2]) + 1],
                columns=self.x[:: int(self.h_vec[-2]) + 1],
            )
        else:
            df = pd.DataFrame(self.z, index=self.y, columns=self.x)

        data = base64.b64encode(df.to_csv().encode()).decode()
        filename = f"INS_structure_factor_{random_number}.csv"

        # model_state for template jinja plot script
        model_state = self.get_model_state()
        model_state["ylabel"] = self.ylabel
        if hasattr(self, "xlabel"):
            model_state["xlabel"] = self.xlabel
        if hasattr(self, "Q0_vec"):
            model_state["Q0"] = self.Q0_vec
        model_state["spectrum_type"] = self.spectrum_type
        if hasattr(self, "ticks_labels"):
            model_state["ticks_positions"] = self.ticks_positions
            model_state["ticks_labels"] = self.ticks_labels
        model_state["filename"] = filename
        model_state["cmap"] = "cividis"

        plotting_script = generate_from_template(model_state)
        plotting_script_data = base64.b64encode(plotting_script.encode()).decode()
        plotting_script_filename = f"plot_script_{random_number}.py"
        return [(data, filename), (plotting_script_data, plotting_script_filename)]

    def _download_data(self, _=None):
        packed_data = self.prepare_data_for_download()
        for data, filename in packed_data:
            self._download(data, filename)

    @staticmethod
    def _download(payload, filename):
        from IPython.display import Javascript

        javas = Javascript(
            """
            var link = document.createElement('a');
            link.href = 'data:text;charset=utf-8;base64,{payload}'
            link.download = "{filename}"
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            """.format(payload=payload, filename=filename)
        )
        display(javas)

    def _clone(self):
        # in case we want to clone the model.
        # This is the case when we have the same data and we inject in three
        # different models: we don't need to fetch three times.
        return copy.deepcopy(self)
