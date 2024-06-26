from aiida.orm import load_code, Dict
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType

from aiida_vibroscopy.common.properties import PhononProperty

HarmonicWorkChain = WorkflowFactory("vibroscopy.phonons.harmonic")

"""
The logic is that HarmonicWorkchain can run PhononWorkchain and DielectricWorkchain, skipping the second
but not the firstfor now we do not support to run only DielectricWorkchain.
"""


def get_builder(codes, structure, parameters):
    from copy import deepcopy

    protocol = parameters["workchain"].pop("protocol", "fast")
    pw_code = codes.get("pw")
    phonopy_code = codes.get("phonopy")
    phonon_property = PhononProperty[
        parameters["harmonic"].pop("phonon_property", "none")
    ]
    supercell_matrix = parameters["harmonic"].pop("supercell_selector", None)
    dielectric_property = parameters["harmonic"].pop("dielectric_property", "none")
    polar = parameters["harmonic"].pop("material_is_polar", "off")

    if polar == "on":
        dielectric_property = "raman"

    scf_overrides = deepcopy(parameters["advanced"])
    overrides = {
        "phonon": {
            "scf": scf_overrides,
            "supercell_matrix": supercell_matrix,
        },
        "dielectric": {"scf": scf_overrides, "property": dielectric_property},
    }

    builder = HarmonicWorkChain.get_builder_from_protocol(
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

    builder.phonon.phonopy.metadata.options.resources = (
        builder.phonopy.metadata.options.resources
    )

    # should be automatic inside HarmonicWorkchain.
    builder.phonon.phonopy.parameters = Dict(dict={})
    builder.phonopy.parameters = builder.phonon.phonopy.parameters
    builder.phonon.phonopy.code = builder.phonopy.code

    builder.phonopy.parameters = Dict(dict=phonon_property.value)

    return builder


def trigger_workchain(name, parameters):
    if name not in ["harmonic", "phonons", "dielectric", "iraman"]:
        return True
    import copy

    parameters_ = copy.deepcopy(parameters)
    harmonic_params = parameters_.pop("harmonic", {})
    phonon_property = harmonic_params.pop("phonon_property", "none")
    dielectric_property = harmonic_params.pop("dielectric_property", "none")
    polar = harmonic_params.pop("material_is_polar", "off")

    trigger_spectrum = parameters_.pop("iraman", {}).pop("spectrum", False)

    trigger_harmonic = (
        polar == "on" or (phonon_property != "none" and dielectric_property != "none")
    ) and not trigger_spectrum
    trigger_phonons = (
        phonon_property != "none"
        and dielectric_property == "none"
        and not trigger_harmonic
    ) and not trigger_spectrum
    trigger_dielectric = (
        phonon_property == "none"
        and dielectric_property != "none"
        and not trigger_harmonic
    ) and not trigger_spectrum

    trigger = {
        "iraman": trigger_spectrum,
        "harmonic": trigger_harmonic,
        "phonons": trigger_phonons,
        "dielectric": trigger_dielectric,
    }

    return trigger[name]


workchain_and_builder = {
    "workchain": HarmonicWorkChain,
    "exclude": ("clean_workdir",),
    "get_builder": get_builder,
}
