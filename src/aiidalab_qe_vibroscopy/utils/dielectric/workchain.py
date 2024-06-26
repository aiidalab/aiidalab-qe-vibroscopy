from aiida.orm import load_code, Str
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType

DielectricWorkChain = WorkflowFactory("vibroscopy.dielectric")

"""
The logic is that HarmonicWorkchain can run PhononWorkchain and DielectricWorkchain, skipping the second
but not the first: so we add also the possibility to run only DielectricWorkchain.
"""


def get_builder(codes, structure, parameters):
    protocol = parameters["workchain"].pop("protocol", "fast")
    pw_code = codes.get("pw")
    dielectric_property = parameters["harmonic"].pop("dielectric_property", "none")

    """
    here we set a readable input anyway, even if we do not run this workflow
    but only the HarmonicWorkchain for phonons-only calculations.
    The problem is that if we select both Dielectric and Phonons calculations,
    this protocol is called, even if never used.
    """
    if dielectric_property == "none":
        dielectric_property = "dielectric"

    overrides = {
        "scf": parameters["advanced"],
    }

    builder = DielectricWorkChain.get_builder_from_protocol(
        code=pw_code,
        structure=structure,
        protocol=protocol,
        overrides=overrides,
        electronic_type=ElectronicType(parameters["workchain"]["electronic_type"]),
        spin_type=SpinType(parameters["workchain"]["spin_type"]),
        initial_magnetic_moments=parameters["advanced"]["initial_magnetic_moments"],
    )

    builder.property = dielectric_property

    builder.pop("clean_workdir", None)

    return builder


workchain_and_builder = {
    "workchain": DielectricWorkChain,
    "exclude": ("clean_workdir",),
    "get_builder": get_builder,
}
