from aiidalab_qe.common.panel import ResultsModel
from aiida.common.extendeddicts import AttributeDict
import traitlets as tl
from aiidalab_qe_vibroscopy.utils.euphonic.data.export_vibronic_to_euphonic import (
    export_euphonic_data,
)


class EuphonicModel(ResultsModel):
    node = tl.Instance(AttributeDict, allow_none=True)

    def fetch_data(self):
        ins_data = export_euphonic_data(self.node)
        self.fc = ins_data["fc"]
        self.q_path = ins_data["q_path"]
