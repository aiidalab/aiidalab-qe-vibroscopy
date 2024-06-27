import pytest

@pytest.mark.usefixtures("sssp")
def test_settings():
    """Test the settings of the vibroscopy app."""
    
    from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep
    from ase.build import bulk
    from aiida.orm import StructureData

    configure_step = ConfigureQeAppWorkChainStep()
    # set the input structure
    silicon = bulk('Si', 'diamond', a=5.43)
    structure = StructureData(ase=silicon)
    configure_step.input_structure = structure
    # select vibrational properties
    configure_step.workchain_settings.properties["vibronic"].run.value = True
    assert configure_step.settings["vibronic"].supercell_widget.layout.display == "block"
    # test get_panel_value for default starting values
    parameters = configure_step.settings["vibronic"].get_panel_value()
    print("parameters", parameters)
    assert parameters["simulation_mode"] == 1
    assert parameters["supercell_selector"] == [2, 2, 2]
    # test the supercell hint for the silicon
    configure_step.settings["vibronic"]._suggest_supercell()
    parameters = configure_step.settings["vibronic"].get_panel_value()
    assert parameters["supercell_selector"] == [4, 4, 4]
    # test that if we select the IR/Raman in unit cell or dielectric properties simulation, we do not have the hint
    configure_step.settings["vibronic"].calc_options.value = 2
    assert configure_step.settings["vibronic"].supercell_widget.layout.display == "none"
    configure_step.settings["vibronic"].calc_options.value = 4
    assert configure_step.settings["vibronic"].supercell_widget.layout.display == "none"
    # test that we go back to display it when we select the other modes (0 and 2)
    configure_step.settings["vibronic"].calc_options.value = 1
    assert configure_step.settings["vibronic"].supercell_widget.layout.display == "block"
    configure_step.settings["vibronic"].calc_options.value = 3
    assert configure_step.settings["vibronic"].supercell_widget.layout.display == "block"
    # test the reset
    configure_step.settings["vibronic"].reset()
    parameters = configure_step.settings["vibronic"].get_panel_value()
    assert parameters["simulation_mode"] == 1
    assert parameters["supercell_selector"] == [2, 2, 2]