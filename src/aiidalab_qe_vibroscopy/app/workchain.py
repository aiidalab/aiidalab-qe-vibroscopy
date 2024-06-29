from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiida_quantumespresso.workflows.pw.bands import PwBaseWorkChain


from aiida import orm

VibroWorkChain = WorkflowFactory("vibroscopy_app.vibro")


def create_resource_config(code_details):
    """
    Create a dictionary with resource configuration based on codes for 'pw'.

    Parameters:
        codes (dict): A dictionary containing the code configuration.

    Returns:
        dict: A nested dictionary with structured resource configurations.
    """
    options =  {
        "options": {
            "resources": {
                "num_machines": code_details["nodes"],
                "num_mpiprocs_per_machine": code_details["ntasks_per_node"],
                "num_cores_per_mpiproc": code_details["cpus_per_task"],
            }, 

        }
    }
    
    if "max_wallclock_seconds" in code_details:
        options["options"]["max_wallclock_seconds"] = code_details["max_wallclock_seconds"]

        
    return options


def get_builder(codes, structure, parameters):
    from copy import deepcopy

    protocol = parameters["workchain"].pop("protocol", "fast")
    pw_code = codes.get("pw")["code"]
    phonopy_code = codes.get("phonopy")["code"]

    simulation_mode = parameters["vibronic"].pop("simulation_mode", 1)

    # Define the supercell matrix
    supercell_matrix = parameters["vibronic"].pop("supercell_selector", None)

    scf_overrides = deepcopy(parameters["advanced"])

    # The following include_all is needed to have forces written
    overrides = {
        "phonon": {
            "scf": scf_overrides,
            "supercell_matrix": supercell_matrix,
        },
        "dielectric": {"scf": scf_overrides},
    }

    # Update code information with resource configurations
    overrides["dielectric"]["scf"]["pw"]["metadata"] = create_resource_config(
        codes.get("pw")
    )
    overrides["phonon"]["scf"]["pw"]["metadata"] = create_resource_config(
        codes.get("pw")
    )

    # Parallelization for phonon calculation
    if "parallelization" in codes.get("pw"):
        overrides["phonon"]["scf"]["pw"]["parallelization"] = orm.Dict(
            codes.get("pw")["parallelization"]
        )

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
        simulation_mode=simulation_mode,
        overrides=overrides,
        electronic_type=ElectronicType(parameters["workchain"]["electronic_type"]),
        spin_type=SpinType(parameters["workchain"]["spin_type"]),
        initial_magnetic_moments=parameters["advanced"]["initial_magnetic_moments"],
    )

    return builder


def update_inputs(inputs, ctx):
    """Update the inputs using context."""
    inputs.structure = ctx.current_structure


workchain_and_builder = {
    "workchain": VibroWorkChain,
    "exclude": ("clean_workdir",),
    "get_builder": get_builder,
    "update_inputs": update_inputs,
}
