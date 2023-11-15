from aiida.orm import load_code, Dict, Bool
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiida_vibroscopy.common.properties import PhononProperty

VibroWorkChain = WorkflowFactory("vibroscopy_app.vibro")

def get_builder(codes, structure, parameters):
    from copy import deepcopy
    protocol = parameters["workchain"].pop("protocol", "fast")
    pw_code = codes.get("pw")
    phonopy_code = codes.get("phonopy")
    
    phonon_property = parameters["vibronic"].pop("phonon_property","none")
    if phonon_property!="none":
        phonon_property = PhononProperty[phonon_property]
    polar = parameters["vibronic"].pop("material_is_polar", "off")
    supercell_matrix = parameters["vibronic"].pop("supercell_selector",None)
        
    dielectric_property = parameters["vibronic"]["dielectric_property"]
    
    spectrum = parameters["vibronic"]["spectrum"]
    
    trigger = "phonon"
    
    if spectrum != 'off':
        trigger = "iraman"
        dielectric_property = spectrum
        
    if polar=="on" and trigger == "phonon":
        #the material is polar, so we need to run HarmonicWChain instead of PhononWChain.
        trigger = "harmonic"
        dielectric_property = "raman"
        
    if trigger not in ["iraman","harmonic"] and (phonon_property != "none" and dielectric_property != "none"):
        trigger="harmonic" 
    if trigger not in ["iraman","harmonic"] and (phonon_property == "none" and dielectric_property != "none"):
        trigger="dielectric"
         
    scf_overrides = deepcopy(parameters["advanced"])
    overrides = {
        "phonon":{
            "scf": scf_overrides,
            "supercell_matrix":supercell_matrix,
            "clean_workdir": Bool(False),
        },
        "dielectric":{
            "scf": scf_overrides,
            "property":dielectric_property,
            "settings.sleep_submission_time": 0.0,
            "clean_workdir": Bool(False),
        },
        "settings":{
        "run_parallel": False
        },
        "clean_workdir": Bool(False),
    }
    
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