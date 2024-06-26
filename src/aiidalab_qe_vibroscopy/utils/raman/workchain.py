from aiida.orm import load_code, Dict
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType

IRamanSpectraWorkChain = WorkflowFactory("vibroscopy.spectra.iraman")


def get_builder(codes, structure, parameters):
    from copy import deepcopy

    protocol = parameters["workchain"].pop("protocol", "fast")
    pw_code = codes.get("pw")
    phonopy_code = codes.get("phonopy")
    supercell_matrix = parameters["iraman"].pop("supercell_selector", None)
    dielectric_property = parameters["iraman"]["spectrum"]
    res = {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 1,
    }

    scf_overrides = deepcopy(parameters["advanced"])
    overrides = {
        "phonon": {
            "scf": scf_overrides,
            "supercell_matrix": supercell_matrix,
        },
        "dielectric": {"scf": scf_overrides, "property": dielectric_property},
    }

    builder = IRamanSpectraWorkChain.get_builder_from_protocol(
        code=pw_code,
        structure=structure,
        protocol=protocol,
        overrides=overrides,
        electronic_type=ElectronicType(parameters["workchain"]["electronic_type"]),
        spin_type=SpinType(parameters["workchain"]["spin_type"]),
        initial_magnetic_moments=parameters["advanced"]["initial_magnetic_moments"],
    )

    builder.dielectric.property = dielectric_property

    builder.phonon.phonopy.code = phonopy_code
    builder.phonon.phonopy.metadata.options.resources = res
    builder.phonon.phonopy.parameters = Dict({})

    return builder


workchain_and_builder = {
    "workchain": IRamanSpectraWorkChain,
    "exclude": ("clean_workdir",),
    "get_builder": get_builder,
}
