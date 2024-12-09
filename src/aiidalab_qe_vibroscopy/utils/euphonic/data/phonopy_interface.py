import tempfile
import pathlib
import base64
from typing import Optional
import euphonic
from phonopy.file_IO import write_force_constants_to_hdf5

from aiidalab_qe_vibroscopy.utils.euphonic.data.structure_factors import (
    blockPrint,
    enablePrint,
)


def generate_force_constant_from_phonopy(
    phonopy_calc=None,
    path: str = None,
    summary_name: str = None,
    born_name: Optional[str] = None,
    fc_name: str = "FORCE_CONSTANTS",
    fc_format: Optional[str] = None,
    mode="stream",  # "download" to have the download of phonopy.yaml and fc.hdf5 . TOBE IMPLEMENTED.
):
    """
    Basically allows to obtain the ForceConstants instance from phonopy, both via files (from the second
    input parameters we have the same one of `euphonic.ForceConstants.from_phonopy`), or via a
    PhonopyCalculation instance. Respectively, the two ways will support independent euphonic app and integration
    of Euphonic into aiidalab.
    """
    blockPrint()

    ####### This is done to support the detached app (from aiidalab) with the same code:
    if path and summary_name:
        fc = euphonic.ForceConstants.from_phonopy(
            path=path,
            summary_name=summary_name,
            fc_name=fc_name,
        )
        return fc
    elif not phonopy_calc:
        raise NotImplementedError(
            "Please provide or the files or the phonopy calculation node."
        )

    ####### This is almost copied from PhonopyCalculation and is done to support functionalities in aiidalab env:
    from phonopy.interface.phonopy_yaml import PhonopyYaml

    kwargs = {}

    if "settings" in phonopy_calc.inputs:
        the_settings = phonopy_calc.inputs.settings.get_dict()
        for key in ["symmetrize_nac", "factor_nac", "subtract_residual_forces"]:
            if key in the_settings:
                kwargs.update({key: the_settings[key]})

    if "phonopy_data" in phonopy_calc.inputs:
        ph = phonopy_calc.inputs.phonopy_data.get_phonopy_instance(**kwargs)
        p2s_map = phonopy_calc.inputs.phonopy_data.get_cells_mappings()["primitive"][
            "p2s_map"
        ]
        ph.produce_force_constants()
    elif "force_constants" in phonopy_calc.inputs:
        ph = phonopy_calc.inputs.force_constants.get_phonopy_instance(**kwargs)
        p2s_map = phonopy_calc.inputs.force_constants.get_cells_mappings()["primitive"][
            "p2s_map"
        ]
        ph.force_constants = phonopy_calc.inputs.force_constants.get_array(
            "force_constants"
        )

    #######

    # Create temporary directory
    #
    with tempfile.TemporaryDirectory() as dirpath:
        # phonopy.yaml generation:
        phpy_yaml = PhonopyYaml()
        phpy_yaml.set_phonon_info(ph)
        phpy_yaml_txt = str(phpy_yaml)

        with open(
            pathlib.Path(dirpath) / "phonopy.yaml", "w", encoding="utf8"
        ) as handle:
            handle.write(phpy_yaml_txt)

        # Force constants hdf5 file generation:
        # all this is needed to load the euphonic instance, in case no FC are written in phonopy.yaml
        # which is the case

        write_force_constants_to_hdf5(
            force_constants=ph.force_constants,
            filename=pathlib.Path(dirpath) / "fc.hdf5",
            p2s_map=p2s_map,
        )

        # Here below we trigger the download mode. Can be improved avoiding the repetitions of lines
        if mode == "download":
            with open(
                pathlib.Path(dirpath) / "phonopy.yaml", "r", encoding="utf8"
            ) as handle:
                file_content = handle.read()
                phonopy_yaml_bitstream = base64.b64encode(file_content.encode()).decode(
                    "utf-8"
                )

            with open(
                pathlib.Path(dirpath) / "fc.hdf5",
                "rb",
            ) as handle:
                file_content = handle.read()
                fc_hdf5_bitstream = base64.b64encode(file_content).decode()

            return phonopy_yaml_bitstream, fc_hdf5_bitstream

        # Read force constants (fc.hdf5) and summary+NAC (phonopy.yaml)

        fc = euphonic.ForceConstants.from_phonopy(
            path=dirpath,
            summary_name="phonopy.yaml",
            fc_name="fc.hdf5",
        )
        # print(filename)
        # print(dirpath)
    enablePrint()
    return fc
