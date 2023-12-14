from aiida.orm import load_code, Dict
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiida_vibroscopy.common.properties import PhononProperty
from aiida_quantumespresso.workflows.pw.bands import PwBaseWorkChain

VibroWorkChain = WorkflowFactory("vibroscopy_app.vibro")


def get_builder(codes, structure, parameters):
    from copy import deepcopy

    protocol = parameters["workchain"].pop("protocol", "fast")
    pw_code = codes.get("pw")
    phonopy_code = codes.get("phonopy")

    calc_option = parameters["vibronic"].pop("calc_options", "raman")

    calc_option_mapping = {
        "ph_bands": PhononProperty.BANDS,
        "ph_pdos": PhononProperty.DOS,
        "ph_therm": PhononProperty.THERMODYNAMIC,
    }

    # Define the phonon property to be calculated
    phonon_property = calc_option_mapping.get(calc_option, PhononProperty.NONE)

    # Define if the material is polar
    polar = parameters["vibronic"].pop("material_is_polar", "off")
    # Define the supercell matrix
    supercell_matrix = parameters["vibronic"].pop("supercell_selector", None)

    # initialize the dielectric property
    dielectric_property = "none"
    # Logic to define the trigger and the dielectric property
    if calc_option in ["ph_bands", "ph_pdos", "ph_therm"]:
        if polar == "on":
            trigger = "harmonic"
            dielectric_property = "raman"
        else:
            trigger = "phonon"
    elif calc_option in ["raman", "ir"]:
        trigger = "iraman"
        dielectric_property = calc_option
    elif calc_option == "dielectric":
        trigger = "dielectric"
        dielectric_property = "dielectric"

    scf_overrides = deepcopy(parameters["advanced"])
    overrides = {
        "phonon": {
            "scf": scf_overrides,
            "supercell_matrix": supercell_matrix,
        },
        "dielectric": {"scf": scf_overrides, "property": dielectric_property},
    }

    # Only for 2D and 1D materials
    if structure.pbc != (True, True, True):
        if "kpoints_distance" not in parameters["advanced"]:
            overrides["dielectric"]["scf"][
                "kpoints_distance"
            ] = PwBaseWorkChain.get_protocol_inputs(protocol)["kpoints_distance"]

    builder = VibroWorkChain.get_builder_from_protocol(
        pw_code=pw_code,
        phonopy_code=phonopy_code,
        structure=structure,
        protocol=protocol,
        trigger=trigger,
        phonon_property=phonon_property,
        dielectric_property=dielectric_property,
        overrides=overrides,
        electronic_type=ElectronicType(parameters["workchain"]["electronic_type"]),
        spin_type=SpinType(parameters["workchain"]["spin_type"]),
        initial_magnetic_moments=parameters["advanced"]["initial_magnetic_moments"],
    )

    return builder


workchain_and_builder = {
    "workchain": VibroWorkChain,
    "exclude": ("clean_workdir",),
    "get_builder": get_builder,
}
