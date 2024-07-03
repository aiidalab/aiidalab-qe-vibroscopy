"""Implementation of the VibroWorkchain for managing the aiida-vibroscopy workchains."""
from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import ToContext, WorkChain, calcfunction
from aiida.orm import AbstractCode, Int, Float, Dict, Code, StructureData, load_code
from aiida.plugins import WorkflowFactory, CalculationFactory
from aiida_quantumespresso.utils.mapping import prepare_process_inputs
from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiida.engine import WorkChain, calcfunction, if_
from aiida_vibroscopy.common.properties import PhononProperty
from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
import math
import numpy as np


GAMMA = "$\Gamma$"

IRamanSpectraWorkChain = WorkflowFactory("vibroscopy.spectra.iraman")
HarmonicWorkChain = WorkflowFactory("vibroscopy.phonons.harmonic")
DielectricWorkChain = WorkflowFactory("vibroscopy.dielectric")
PhononWorkChain = WorkflowFactory("vibroscopy.phonons.phonon")
PhonopyCalculation = CalculationFactory("phonopy.phonopy")


def generate_2d_path(symmetry_type, eta=None, nu=None):
    PATH_SYMMETRY_2D = {
        "hexagonal": {
            "band": [0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.3333, 0.3333, 0.0, 0.0, 0.0, 0.0],
            "labels": [GAMMA, "$\\mathrm{M}$", "$\\mathrm{K}$", GAMMA],
        },
        "square": {
            "band": [0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0],
            "labels": [GAMMA, "$\\mathrm{X}$", "$\\mathrm{M}$", GAMMA],
        },
        "rectangular": {
            "band": [
                0.0,
                0.0,
                0.0,
                0.5,
                0.0,
                0.0,
                0.5,
                0.5,
                0.0,
                0.0,
                0.5,
                0.0,
                1.0,
                0.0,
                0.0,
            ],
            "labels": [GAMMA, "$\\mathrm{X}$", "$\\mathrm{S}$", "$\\mathrm{Y}$", GAMMA],
        },
        "rectangular_centered": {
            "band": [
                0.0,
                0.0,
                0.0,
                0.5,
                0.0,
                0.0,
                1 - eta,
                nu,
                0,
                0.5,
                0.5,
                0.0,
                eta,
                1 - nu,
                0.0,
                1.0,
                0.0,
                0.0,
            ],
            "labels": [
                GAMMA,
                "$\\mathrm{X}$",
                "$\\mathrm{H_1}$",
                "$\\mathrm{C}$",
                "$\\mathrm{H}$",
                GAMMA,
            ],
        },
        "oblique": {
            "band": [
                0.0,
                0.0,
                0.0,
                0.5,
                0.0,
                0.0,
                1 - eta,
                nu,
                0,
                0.5,
                0.5,
                0.0,
                eta,
                1 - nu,
                0.0,
                0.0,
                0.5,
                0.0,
                1.0,
                0.0,
                0.0,
            ],
            "labels": [
                GAMMA,
                "$\\mathrm{X}$",
                "$\\mathrm{H_1}$",
                "$\\mathrm{C}$",
                "$\\mathrm{H}$",
                "$\\mathrm{Y}$",
                GAMMA,
            ],
        },
    }

    if symmetry_type in PATH_SYMMETRY_2D:
        return PATH_SYMMETRY_2D[symmetry_type]
    else:
        raise ValueError("Invalid symmetry type")


def determine_symmetry_path(structure):
    # Tolerance for checking equality
    cell_lengths = structure.cell_lengths
    cell_angles = structure.cell_angles
    tolerance = 1e-3

    # Define symmetry conditions and their corresponding types in a dictionary
    symmetry_conditions = {
        (
            math.isclose(cell_lengths[0], cell_lengths[1], abs_tol=tolerance)
            and math.isclose(cell_angles[2], 120.0, abs_tol=tolerance)
        ): "hexagonal",
        (
            math.isclose(cell_lengths[0], cell_lengths[1], abs_tol=tolerance)
            and math.isclose(cell_angles[2], 90.0, abs_tol=tolerance)
        ): "square",
        (
            math.isclose(cell_lengths[0], cell_lengths[1], abs_tol=tolerance) == False
            and math.isclose(cell_angles[2], 90.0, abs_tol=tolerance)
        ): "rectangular",
        (
            math.isclose(
                cell_lengths[1] * math.cos(math.radians(cell_angles[2])),
                cell_lengths[0] / 2,
                abs_tol=tolerance,
            )
        ): "rectangular_centered",
        (
            math.isclose(cell_lengths[0], cell_lengths[1], abs_tol=tolerance) == False
            and math.isclose(cell_angles[2], 90.0, abs_tol=tolerance) == False
        ): "oblique",
    }

    # Check for symmetry type based on conditions
    for condition, symmetry_type in symmetry_conditions.items():
        if condition:
            if symmetry_type == "rectangular_centered" or "oblique":
                cos_gamma = np.array(structure.cell[0]).dot(structure.cell[1]) / (
                    cell_lengths[0] * cell_lengths[1]
                )
                gamma = np.arccos(cos_gamma)
                eta = (1 - (cell_lengths[0] / cell_lengths[1]) * cos_gamma) / (
                    2 * np.power(np.sin(gamma), 2)
                )
                nu = 0.5 - (eta * cell_lengths[1] * cos_gamma) / cell_lengths[0]
                return generate_2d_path(symmetry_type, eta, nu)

            return generate_2d_path(symmetry_type)
        else:
            raise ValueError("Invalid symmetry type")


class VibroWorkChain(WorkChain):
    "WorkChain to compute vibrational property of a crystal."
    label = "vibro"

    @classmethod
    def define(cls, spec):
        """Specify inputs and outputs."""
        super().define(spec)
        spec.input(
            "structure", valid_type=StructureData
        )  # Maybe not needed as input... just in the protocols. but in this way it is not easy to automate it in the app, after the relaxation. So let's keep it for now.

        spec.expose_inputs(
            PhononWorkChain,
            namespace="phonon",
            exclude=("clean_workdir"),  # AAA check this... maybe not needed.
            namespace_options={
                "required": False,
                "populate_defaults": False,
                "help": 'Inputs for the `PhononWorkChain`, triggered selecting "trigger"="phonon".',
            },
        )
        spec.expose_inputs(
            DielectricWorkChain,
            namespace="dielectric",
            exclude=("clean_workdir"),  # AAA check this... maybe not needed.
            namespace_options={
                "required": False,
                "populate_defaults": False,
                "help": (
                    "Inputs for the `DielectricWorkChain` that will be"
                    "used to calculate the mixed derivatives with electric field."
                ),
            },
            # exclude=('symmetry')
        )
        spec.expose_inputs(
            HarmonicWorkChain,
            namespace="harmonic",
            exclude=("clean_workdir"),  # AAA check this... maybe not needed.
            namespace_options={
                "required": False,
                "populate_defaults": False,
                "help": "Inputs for the `HarmonicWorkChain`.",
            },
        )
        spec.expose_inputs(
            IRamanSpectraWorkChain,
            namespace="iraman",
            exclude=("clean_workdir"),  # AAA check this... maybe not needed.
            namespace_options={
                "required": False,
                "populate_defaults": False,
                "help": "Inputs for the `IRamanSpectraWorkChain`.",
            },
        )
        spec.expose_inputs(
            PhonopyCalculation, namespace='phonopy_calc',
            namespace_options={
                'required': False, 'populate_defaults': False,
                'help': (
                    'Inputs for the `PhonopyCalculation` that will'
                    'be used to calculate the inter-atomic force constants, or for post-processing.'
                )
            },
            exclude=['phonopy_data', 'force_constants', 'parameters'],
        )
        spec.input(
            "phonopy_bands_dict",
            valid_type=Dict,
            required=False,
            help="Settings for phonopy bands calculation.",
        )
        spec.input(
            "phonopy_pdos_dict",
            valid_type=Dict,
            required=False,
            help="Settings for phonopy pdos calculation.",
        )
        spec.input(
            "phonopy_thermo_dict",
            valid_type=Dict,
            required=False,
            help="Settings for phonopy pdos calculation.",
        )
        ###
        spec.outline(
            cls.setup,
            cls.vibrate,
            if_(cls.should_run_phonopy)(
                cls.run_phonopy,
            ),
            cls.results,
        )
        ###
        spec.expose_outputs(
            PhononWorkChain,
            exclude=("supercells", "supercells_forces"),
            namespace="phonon",
            namespace_options={
                "required": False,
                "help": "Outputs of the `PhononWorkChain`.",
            },
        )
        spec.expose_outputs(
            DielectricWorkChain,
            exclude=("fields_data"),
            namespace="dielectric",
            namespace_options={
                "required": False,
                "help": "Outputs of the `DielectricWorkChain`.",
            },
        )
        spec.expose_outputs(
            HarmonicWorkChain,
            exclude=("output_phonon", "output_dielectric"),
            namespace="harmonic",
            namespace_options={
                "required": False,
                "help": "Outputs of the `HarmonicWorkChain`.",
            },
        )
        spec.expose_outputs(
            IRamanSpectraWorkChain,
            exclude=("output_phonon", "output_dielectric"),
            namespace="iraman",
            namespace_options={
                "required": False,
                "help": "Outputs of the `IRamanSpectraWorkChain`.",
            },
        )
        spec.output(
            "phonon_bands",
            valid_type=orm.BandsData,
            required=False,
            help="Calculated phonon band structure.",
        )
        spec.output(
            "phonon_pdos",
            valid_type=orm.XyData,
            required=False,
            help="Calculated projected DOS.",
        )
        spec.output(
            "phonon_thermo",
            valid_type=orm.XyData,
            required=False,
            help="Calculated thermal properties.",
        )
        ###
        spec.exit_code(400, "ERROR_WORKCHAIN_FAILED", message="The workchain failed.")

    @classmethod
    def get_builder_from_protocol(
        cls,
        pw_code,
        structure,
        protocol=None,
        phonopy_code=None,
        overrides=None,
        options=None,
        simulation_mode=None,
        phonon_property=PhononProperty.THERMODYNAMIC,
        dielectric_property="raman",
        **kwargs,
    ):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param pw_code: the ``Code`` instance configured for the ``quantumespresso.pw`` plugin.
        :param structure: the ``StructureData`` instance to use.
        :param phonopy_code: the ``Code`` instance configured for the ``phonopy.phonopy`` plugin.
        :param protocol: protocol to use, if not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the protocol.
        :param simulation_mode: what type of simulation to run. Refer to the settings.py of the app.
        :param options: A dictionary of options that will be recursively set for the ``metadata.options`` input for the builder of the pw workchains.
        :param kwargs: additional keyword arguments that will be passed to the ``get_builder_from_protocol`` of all the
            sub processes that are called by this workchain.
        :return: a process builder instance with all inputs defined ready for launch.
        """
        from aiida_quantumespresso.workflows.protocols.utils import recursive_merge

        if simulation_mode not in range(1, 5):
            raise ValueError("simulation_mode not in [1,2,3,4]")

        builder = cls.get_builder()

        if simulation_mode == 1:
            # Running the full workchain: IR/Raman, Phonons, INS, Dielectric

            builder_harmonic = HarmonicWorkChain.get_builder_from_protocol(
                pw_code=pw_code,
                phonopy_code=phonopy_code,
                structure=structure,
                protocol=protocol,
                overrides=overrides,
                phonon_property=phonon_property,
                **kwargs,
            )

            # MB guesses phonopy will always run serially, otherwise choose phono3py
            # also this is needed to be set here.
            builder_harmonic.phonopy.metadata.options.resources = {
                "num_machines": 1,
                "num_mpiprocs_per_machine": 1,
            }

            builder_harmonic.phonon.phonopy.metadata.options.resources = (
                builder_harmonic.phonopy.metadata.options.resources
            )

            # should be automatic inside HarmonicWorkchain.
            builder_harmonic.phonon.phonopy.parameters = Dict(dict={})
            builder_harmonic.phonopy.parameters = (
                builder_harmonic.phonon.phonopy.parameters
            )
            builder_harmonic.phonon.phonopy.code = builder_harmonic.phonopy.code

            builder_harmonic.phonopy.parameters = Dict(dict=phonon_property.value)

            # Setting the `raman` dielectric property, to compute up to the third order derivative wrt finite electric fields.
            builder_harmonic.dielectric.property = dielectric_property

            # To run euphonic: we should be able to get rid of this, using phonopy API.
            builder_harmonic.phonon.phonopy.settings = Dict(
                dict={"keep_phonopy_yaml": True}
            )

            builder.harmonic = builder_harmonic

            # Adding the bands and pdos inputs.
            if structure.pbc != (True, True, True):
                # Generate Mesh for 1D and 2D materials
                builder.harmonic.dielectric.pop("kpoints_parallel_distance", None)
                inputs = {
                    "structure": structure,
                    "distance": orm.Float(0.01),
                    "force_parity": orm.Bool(False),
                    "metadata": {"call_link_label": "create_kpoints_from_distance"},
                }
                kpoints = create_kpoints_from_distance(**inputs)

                builder.phonopy_pdos_dict = Dict(
                    dict={
                        "pdos": "auto",
                        "mesh": kpoints.get_kpoints_mesh()[0],
                        "write_mesh": False,
                    }
                )
                builder.phonopy_thermo_dict = Dict(dict=PhononProperty.THERMODYNAMIC.value)

                if structure.pbc == (True, False, False):
                    builder.phonopy_bands_dict = Dict(
                        dict={
                            "symmetry_tolerance": 1e-3,
                            "band": [0, 0, 0, 1 / 2, 0, 0],
                            "band_points": 100,
                            "band_labels": [GAMMA, "$\\mathrm{X}$"],
                        }
                    )
                    #change symprec for 1D materials to 1e-3
                    builder.phonopy_pdos_dict["symmetry_tolerance"] = 1e-3
                    builder.phonopy_thermo_dict["symmetry_tolerance"] = 1e-3
                    builder.harmonic.symmetry.symprec = orm.Float(1e-3)
                    builder.harmonic.phonon.phonopy.parameters = orm.Dict({"symmetry_tolerance": 1e-3})

                elif structure.pbc == (True, True, False):
                    symmetry_path = determine_symmetry_path(structure)
                    builder.phonopy_bands_dict = Dict(
                        dict={
                            "band": symmetry_path["band"],
                            "band_points": 100,
                            "band_labels": symmetry_path["labels"],
                            "primitive_axes": [
                                1.0,
                                0.0,
                                0.0,
                                0.0,
                                1.0,
                                0.0,
                                0.0,
                                0.0,
                                1.0,
                            ],
                        }
                    )
            else:
                builder.phonopy_bands_dict = Dict(dict=PhononProperty.BANDS.value)
                builder.phonopy_pdos_dict = Dict(
                    dict={
                        "pdos": "auto",
                        "mesh": 150,  # 1000 is too heavy
                        "write_mesh": False,
                    }
                )

                builder.phonopy_thermo_dict = Dict(dict=PhononProperty.THERMODYNAMIC.value)

        elif simulation_mode == 2:

            builder_iraman = IRamanSpectraWorkChain.get_builder_from_protocol(
                code=pw_code,
                structure=structure,
                protocol=protocol,
                overrides=overrides,
                **kwargs,
            )

            builder_iraman.dielectric.property = dielectric_property

            builder_iraman.phonon.phonopy.code = phonopy_code
            builder_iraman.phonon.phonopy.metadata.options.resources = {
                "num_machines": 1,
                "num_mpiprocs_per_machine": 1,
            }
            builder_iraman.phonon.phonopy.parameters = Dict({})

            # Remove kpoints_parallel_distance if for 1D and 2D materials
            if structure.pbc != (True, True, True):
                builder_iraman.dielectric.pop("kpoints_parallel_distance", None)

                if structure.pbc == (True, False, False):
                    builder_iraman.symmetry.symprec = orm.Float(1e-3)
                    builder_iraman.phonon.phonopy.parameters = orm.Dict({"symmetry_tolerance": 1e-3})

            builder.iraman = builder_iraman

        elif simulation_mode == 3:

            builder_phonon = PhononWorkChain.get_builder_from_protocol(
                pw_code=pw_code,
                phonopy_code=phonopy_code,
                structure=structure,
                protocol=protocol,
                overrides=overrides["phonon"],
                phonon_property=phonon_property,
                **kwargs,
            )

            builder_phonon.phonopy.metadata.options.resources = {
                "num_machines": 1,
                "num_mpiprocs_per_machine": 1,
            }

            # To run euphonic: we should be able to get rid of this, using phonopy API.
            builder_phonon.phonopy.settings = Dict(dict={"keep_phonopy_yaml": True})

            # MBO: I do not understand why I have to do this, but it works
            symmetry = builder_phonon.pop("symmetry")
            builder.phonon = builder_phonon
            builder.phonon.symmetry = symmetry

            # Adding the bands and pdos inputs.
            if structure.pbc != (True, True, True):
                # Generate Mesh for 1D and 2D materials
                inputs = {
                    "structure": structure,
                    "distance": orm.Float(0.01),
                    "force_parity": orm.Bool(False),
                    "metadata": {"call_link_label": "create_kpoints_from_distance"},
                }
                kpoints = create_kpoints_from_distance(**inputs)

                builder.phonopy_pdos_dict = Dict(
                    dict={
                        "pdos": "auto",
                        "mesh": kpoints.get_kpoints_mesh()[0],
                        "write_mesh": False,
                    }
                )

                builder.phonopy_thermo_dict = Dict(dict=PhononProperty.THERMODYNAMIC.value)
                
                if structure.pbc == (True, False, False):
                    builder.phonopy_bands_dict = Dict(
                        dict={
                            "symmetry_tolerance": 1e-3,
                            "band": [0, 0, 0, 1 / 2, 0, 0],
                            "band_points": 100,
                            "band_labels": [GAMMA, "$\\mathrm{X}$"],
                        }
                    )
                    #change symprec for 1D materials to 1e-3
                    builder.phonopy_pdos_dict["symmetry_tolerance"] = 1e-3
                    builder.phonopy_thermo_dict["symmetry_tolerance"] = 1e-3
                    builder.phonon.symmetry.symprec = orm.Float(1e-3)
                    builder.phonon.phonopy.parameters = orm.Dict({"symmetry_tolerance": 1e-3})

                elif structure.pbc == (True, True, False):
                    symmetry_path = determine_symmetry_path(structure)
                    builder.phonopy_bands_dict = Dict(
                        dict={
                            "band": symmetry_path["band"],
                            "band_points": 100,
                            "band_labels": symmetry_path["labels"],
                            "primitive_axes": [
                                1.0,
                                0.0,
                                0.0,
                                0.0,
                                1.0,
                                0.0,
                                0.0,
                                0.0,
                                1.0,
                            ],
                        }
                    )
            else:
                builder.phonopy_bands_dict = Dict(dict=PhononProperty.BANDS.value)
                builder.phonopy_pdos_dict = Dict(
                    dict={
                        "pdos": "auto",
                        "mesh": 150,  # 1000 is too heavy
                        "write_mesh": False,
                    }
                )

                builder.phonopy_thermo_dict = Dict(dict=PhononProperty.THERMODYNAMIC.value)

        elif simulation_mode == 4:

            builder_dielectric = DielectricWorkChain.get_builder_from_protocol(
                code=pw_code,
                structure=structure,
                protocol=protocol,
                overrides=overrides["dielectric"],
                **kwargs,
            )

            # MBO: I do not understand why I have to do this, but it works. maybe related with excludes.
            symmetry = builder_dielectric.pop("symmetry")
            if structure.pbc != (True, True, True):
                builder_dielectric.dielectric.pop("kpoints_parallel_distance", None)
            builder.dielectric = builder_dielectric
            builder.dielectric.symmetry = symmetry
            builder.dielectric.property = dielectric_property

        # Deleting the not needed parts of the builder:
        available_wchains = [
            "harmonic",
            "iraman",
            "phonon",
            "dielectric",
        ]
        for wchain_idx in range(1, 5):
            if simulation_mode != wchain_idx:
                builder.pop(available_wchains[wchain_idx - 1], None)
        builder.phonopy_calc.code = phonopy_code
        builder.phonopy_calc.metadata.options.resources = {
                "num_machines": 1,
                "num_mpiprocs_per_machine": 1,
            }
        if structure.pbc == (True, False, False):
            builder.phonopy_calc.parameters = orm.Dict({"symmetry_tolerance": 1e-3})
        builder.structure = structure

        return builder

    def setup(self):
        # setup general contest variables... see in HarmonicWorkChain.
        # also see self.ctx.key = 'phonon'... maybe you can initialise here to simplify next functions.
        # key, class, outputs namespace.
        if "phonon" in self.inputs:
            self.ctx.key = "phonon"
            self.ctx.workchain = PhononWorkChain
            self.ctx.phonopy = (
                self.inputs.phonon.phonopy
            )  # in the ctx, because then I will delete them in the workchain run
        elif "dielectric" in self.inputs:
            self.ctx.key = "dielectric"
            self.ctx.workchain = DielectricWorkChain
        elif "harmonic" in self.inputs:
            self.ctx.key = "harmonic"
            self.ctx.workchain = HarmonicWorkChain
            self.ctx.phonopy = (
                self.inputs.harmonic.phonon.phonopy
            )  # in the ctx, because then I will delete them in the workchain run
        elif "iraman" in self.inputs:
            self.ctx.key = "iraman"
            self.ctx.workchain = IRamanSpectraWorkChain

        self.ctx.structure = self.inputs.structure

    def vibrate(self):
        """Run a WorkChain for vibrational properties."""
        # maybe we can unify this, thanks to a wise setup.
        inputs = AttributeDict(
            self.exposed_inputs(self.ctx.workchain, namespace=self.ctx.key)
        )
        inputs.metadata.call_link_label = self.ctx.key
        inputs.pop("phonopy", None)  # I will run phonopy later in the outline.
        if self.ctx.key in ["phonon", "dielectric"]:
            inputs.scf.pw.structure = self.ctx.structure
        else:
            inputs.structure = self.ctx.structure

        future = self.submit(self.ctx.workchain, **inputs)
        self.report(f"submitting `WorkChain` <PK={future.pk}>")
        self.to_context(**{self.ctx.key: future})

    def should_run_phonopy(self):
        # Final phonopy is needed for modes 1 and 3;
        # namely, we compute bands, pdos and thermo.
        return "phonopy_bands_dict" in self.inputs

    def run_phonopy(self):
        """Run three `PhonopyCalculation` to get (after the calculations of force constants) bands, pdos, thermodynamic quantities."""

        for calc_type in ["bands", "pdos", "thermo"]:
            key = f"phonopy_{calc_type}_calculation"

            # The following copy is done in order to be able to update the AttributesFrozenDict self.ctx.phonopy
            # try in a verdi shell: from plumpy.utils import AttributesFrozendict
            if self.ctx.key == "phonon":
                # See how this works in PhononWorkChain
                inputs = AttributeDict(self.exposed_inputs(PhonopyCalculation, namespace="phonopy_calc"))
                inputs.phonopy_data = self.ctx[self.ctx.key].outputs.phonopy_data
                inputs.parameters = self.inputs[f"phonopy_{calc_type}_dict"]
                
            elif self.ctx.key == "harmonic":
                # See how this works in HarmonicWorkChain
                inputs = AttributeDict(self.exposed_inputs(PhonopyCalculation, namespace="phonopy_calc"))
                inputs.force_constants  = list(self.ctx[self.ctx.key].outputs["vibrational_data"].values())[-1]
                inputs.parameters = self.inputs[f"phonopy_{calc_type}_dict"]
                
            inputs.metadata.call_link_label = key
            future = self.submit(PhonopyCalculation, **inputs)
            self.report(
                f"submitting `PhonopyCalculation` for {calc_type}, <PK={future.pk}>"
            )
            self.to_context(**{calc_type: future})

    def results(self):
        """Inspect all sub-processes."""
        workchain = self.ctx[self.ctx.key]
        failed = False

        if not workchain.is_finished_ok:
            self.report(f"the child WorkChain with <PK={workchain.pk}> failed")
            failed = True
        else:
            self.out_many(
                self.exposed_outputs(
                    self.ctx[self.ctx.key], self.ctx.workchain, namespace=self.ctx.key
                )
            )

        if "bands" in self.ctx.keys():
            if self.ctx["bands"].is_finished_ok:
                self.out("phonon_bands", self.ctx["bands"].outputs.phonon_bands)
            else:
                self.report(f"the child bands PhonopyCalculation failed")
                failed = True

            if self.ctx["pdos"].is_finished_ok:
                self.out("phonon_pdos", self.ctx["pdos"].outputs.projected_phonon_dos)
            else:
                self.report(f"the child pdos PhonopyCalculation failed")
                failed = True

            if self.ctx["thermo"].is_finished_ok:
                self.out("phonon_thermo", self.ctx["thermo"].outputs.thermal_properties)
            else:
                self.report(f"the child thermo PhonopyCalculation failed")
                failed = True

        if failed:
            return self.exit_codes.ERROR_WORKCHAIN_FAILED
