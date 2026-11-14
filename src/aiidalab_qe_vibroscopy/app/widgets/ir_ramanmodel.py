from __future__ import annotations
from aiidalab_qe.common.mvc import Model
from aiida.common.extendeddicts import AttributeDict
from ase.atoms import Atoms
import traitlets as tl

from aiidalab_qe_vibroscopy.utils.raman.result import export_iramanworkchain_data


class IRRamanModel(Model):
    vibro = tl.Instance(AttributeDict, allow_none=True)
    input_structure = tl.Instance(Atoms, allow_none=True)

    needs_raman_tab = tl.Bool()
    needs_ir_tab = tl.Bool()

    def fetch_data(self):
        spectra_data = export_iramanworkchain_data(self.vibro)
        if spectra_data["Ir"] == "No IR modes detected.":
            self.needs_ir_tab = False
        else:
            self.needs_ir_tab = True
        if spectra_data["Raman"] == "No Raman modes detected.":
            self.needs_raman_tab = False
        else:
            self.needs_raman_tab = True
