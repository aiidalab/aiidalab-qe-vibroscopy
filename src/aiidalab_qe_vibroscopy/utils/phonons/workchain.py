from aiida.orm import load_code
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiida_vibroscopy.common.properties import PhononProperty

PhononWorkChain = WorkflowFactory("vibroscopy.phonons.phonon")


def get_builder(codes, structure, parameters):
    from copy import deepcopy

    pw_code = codes.get("pw")
    phonopy_code = codes.get("phonopy", None)
    phonon_property = PhononProperty[
        parameters["harmonic"].pop("phonon_property", "none")
    ]
    supercell_matrix = parameters["harmonic"].pop("supercell_selector", None)
    protocol = parameters["workchain"].pop("protocol", "fast")
    scf_overrides = deepcopy(parameters["advanced"])
    overrides = {
        "scf": scf_overrides,
        "supercell_matrix": supercell_matrix,
    }

    builder = PhononWorkChain.get_builder_from_protocol(
        pw_code=pw_code,
        phonopy_code=phonopy_code,
        structure=structure,
        protocol=protocol,
        overrides=overrides,
        phonon_property=phonon_property,
        electronic_type=ElectronicType(parameters["workchain"]["electronic_type"]),
        spin_type=SpinType(parameters["workchain"]["spin_type"]),
        initial_magnetic_moments=parameters["advanced"]["initial_magnetic_moments"],
    )

    # MB supposes phonopy will always run serially, otherwise choose phono3py
    # also this is needed to be set here.
    builder.phonopy.metadata.options.resources = {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 1,
    }

    return builder


workchain_and_builder = {
    "workchain": PhononWorkChain,
    "exclude": ("clean_workdir",),
    "get_builder": get_builder,
}
