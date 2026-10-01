from aiidalab_qe.common.mvc import Model
from aiida import orm
from ase.atoms import Atoms
import traitlets as tl


class IRRamanModel(Model):
    vibrational_data_uuid = tl.Unicode()
    input_structure = tl.Instance(Atoms, allow_none=True)
    needs_raman_tab = tl.Bool()
    needs_ir_tab = tl.Bool()

    def fetch_data(self):
        data = orm.load_node(self.vibrational_data_uuid)
        arrays = data.get_arraynames()
        self.needs_ir_tab = (
            "born_charges" in arrays and len(data.run_powder_ir_intensities()[0]) > 0
        )
        self.needs_raman_tab = (
            "raman_tensors" in arrays
            and len(data.run_powder_raman_intensities()[0]) > 0
        )
