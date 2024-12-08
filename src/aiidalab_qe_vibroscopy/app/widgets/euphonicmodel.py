import numpy as np
import traitlets as tl
import copy

from aiidalab_qe_vibroscopy.utils.euphonic.data_manipulation.intensity_maps import (
    AttrDict,
    produce_bands_weigthed_data,
    produce_powder_data,
    generated_curated_data,
    par_dict,
    par_dict_powder,
    export_euphonic_data,
    generate_force_constant_instance,
)

from aiidalab_qe_vibroscopy.utils.euphonic.tab_widgets.euphonic_q_planes_widgets import (
    produce_Q_section_modes,
    produce_Q_section_spectrum,
)

from aiidalab_qe.common.mvc import Model


class EuphonicBaseResultsModel(Model):
    """Model for the neutron scattering results panel."""

    # Here we mode all the model and data-controller, i.e. all the data and their
    # manipulation to produce new spectra.
    # plot-controller, i.e. scales, colors and so on, should be attached to the widget, I think.
    # the above last point should be discussed more.

    # For now, we do only the first of the following:
    # 1. single crystal data: sc
    # 2. powder average: pa
    # 3. Q planes: qp

    spectra = {}
    path = []
    q_path = None
    spectrum_type = "single_crystal"
    x_label = None
    y_label = "Energy (meV)"
    detached_app = False

    # Settings for single crystal and powder average
    q_spacing = tl.Float(0.01)
    energy_broadening = tl.Float(0.05)
    energy_bins = tl.Int(200)
    temperature = tl.Float(0)
    weighting = tl.Unicode("coherent")
    
    def set_model_state(self, parameters: dict):
        for k, v in parameters.items():
            setattr(self, k, v)

    def _get_default(self, trait):
        if trait in ["h_vec", "k_vec"]:
            return [1, 1, 1, 100, 1]
        elif trait == "Q0_vec":
            return [0.0, 0.0, 0.0]
        return self.traits()[trait].default_value

    def get_model_state(self):
        return {trait: getattr(self, trait) for trait in self.traits()}

    def reset(
        self,
    ):
        with self.hold_trait_notifications():
            for trait in self.traits():
                setattr(self, trait, self._get_default(trait))

    def fetch_data(self):
        """Fetch the data from the database or from the uploaded files."""
        # 1. from aiida, so we have the node
        if self.node:
            ins_data = export_euphonic_data(self.node)
            self.fc = ins_data["fc"]
            self.q_path = ins_data["q_path"]
        # 2. from uploaded files...
        else:
            self.fc = self.upload_widget._read_phonopy_files(
                fname=self.fname,
                phonopy_yaml_content=self._model.phonopy_yaml_content,
                fc_hdf5_content=self._model.fc_hdf5_content,
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
        self.add_traits(npts=tl.Int(100))

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
        # and the
        if not hasattr(self, "parameters"):
            self._inject_single_crystal_settings()
            
        if self.spectrum_type == "q_planes":
            self._get_qsection_spectra()
            return

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
        if spectrum_type == "single_crystal":  # single crystal case
            self.x, self.y = np.meshgrid(
                spectra[0].x_data.magnitude, spectra[0].y_data.magnitude
            )
            (
                final_xspectra,
                final_zspectra,
                self.ticks_positions,
                self.ticks_labels,
            ) = generated_curated_data(spectra)
            
            self.z = final_zspectra.T
            self.y = self.y[:,0]
            self.x = None # we have the ticks positions and labels
            
            self.xlabel = ""
            self.ylabel = "Energy (meV)"
        
        elif spectrum_type == "powder":  # powder case
            # Spectrum2D as output of the powder data
            self.x, self.y = np.meshgrid(
                spectra.x_data.magnitude, spectra.y_data.magnitude
            )

            # we don't need to curate the powder data, at variance with the single crystal case.
            # We can directly use them:
            self.x = spectra.x_data.magnitude[0]
            self.y = self.y[:,0]
            self.z = spectra.z_data.magnitude.T

    def _get_qsection_spectra(
        self,
    ):
        # This is used to update the spectra in the case we plot the Q planes (the third tab).
        parameters_qplanes = AttrDict(
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
            h=parameters_qplanes.h,
            k=parameters_qplanes.k,
            Q0=parameters_qplanes.Q0,
            n_h=parameters_qplanes.n_h,
            n_k=parameters_qplanes.n_k,
            h_extension=parameters_qplanes.h_extension,
            k_extension=parameters_qplanes.k_extension,
            temperature=parameters_qplanes.temperature,
        )

        self.av_spec, self.z, self.x, self.y, self.labels = (
            produce_Q_section_spectrum(
                modes,
                q_array,
                h_array,
                k_array,
                ecenter=parameters_qplanes.ecenter,
                deltaE=parameters_qplanes.deltaE,
                bins=parameters_qplanes.bins,
                spectrum_type=parameters_qplanes.spectrum_type,
                dw=dw,
                labels=labels,
            )
        )

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
        phonopy_yaml, fc_hdf5 = generate_force_constant_instance(
            self.node.phonon_bands.creator, mode="download"
        )
        return phonopy_yaml, fc_hdf5
    
    def prepare_data_for_download(self):
        raise NotImplementedError("Need to implement the download of a CSV file")
    
    def _clone(self):
        # in case we want to clone the model. 
        # This is the case when we have the same data and we inject in three
        # different models: we don't need to fetch three times.
        return copy.deepcopy(self)