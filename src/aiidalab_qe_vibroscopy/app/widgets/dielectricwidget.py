import ipywidgets as ipw
from aiidalab_qe_vibroscopy.app.widgets.dielectricmodel import DielectricModel
from aiidalab_widgets_base import LoadingWidget

from aiidalab_qe.common.infobox import InAppGuide


class DielectricWidget(ipw.VBox):
    """
    Widget for displaying dielectric properties results
    """

    def __init__(self, model: DielectricModel, node: None, **kwargs):
        super().__init__(
            children=[LoadingWidget("Loading widgets")],
            **kwargs,
        )
        self._model = model

        self.rendered = False
        self._model.vibro = node

    def render(self):
        if self.rendered:
            return

        self.dielectric_results_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            The DielectricWorkchain computes different properties: <br>
                <em style="display: inline-block; margin-left: 20px;">-High freq. dielectric tensor </em> <br>
                <em style="display: inline-block; margin-left: 20px;">-Born charges </em> <br>
                <em style="display: inline-block; margin-left: 20px;">-Raman tensors </em> <br>
                <em style="display: inline-block; margin-left: 20px;">-The non-linear optical susceptibility tensor </em> <br>
                All information can be downloaded as a JSON file. <br>

            </div>"""
        )

        self.site_selector = ipw.Dropdown(
            layout=ipw.Layout(width="450px"),
            description="Select atom site:",
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self._model, "site_selector_options"),
            (self.site_selector, "options"),
        )
        self.site_selector.observe(self._on_site_change, names="value")

        self.download_button = ipw.Button(
            description="Download Data", icon="download", button_style="primary"
        )

        self.download_button.on_click(self._model.download_data)

        # HTML table with the dielectric tensor
        self.dielectric_tensor_table = ipw.HTML()
        ipw.link(
            (self._model, "dielectric_tensor_table"),
            (self.dielectric_tensor_table, "value"),
        )

        # HTML table with the Born charges @ site
        self.born_charges_table = ipw.HTML()
        ipw.link(
            (self._model, "born_charges_table"),
            (self.born_charges_table, "value"),
        )

        # HTML table with the Raman tensors @ site
        self.raman_tensors_table = ipw.HTML()
        ipw.link(
            (self._model, "raman_tensors_table"),
            (self.raman_tensors_table, "value"),
        )

        self.children = [
            InAppGuide(identifier="dielectric-results"),
            self.dielectric_results_help,
            ipw.HTML("<h3>Dielectric tensor</h3>"),
            self.dielectric_tensor_table,
            self.site_selector,
            ipw.HBox(
                [
                    ipw.VBox(
                        [
                            ipw.HTML("<h3>Born effective charges</h3>"),
                            self.born_charges_table,
                        ]
                    ),
                    ipw.VBox(
                        [
                            ipw.HTML("<h3>Raman tensor </h3>"),
                            self.raman_tensors_table,
                        ]
                    ),
                ]
            ),
            self.download_button,
        ]

        self.rendered = True
        self._initial_view()

    def _initial_view(self):
        self._model.fetch_data()
        self._model.set_initial_values()
        self.dielectric_tensor_table.layout = ipw.Layout(width="300px", height="auto")
        self.born_charges_table.layout = ipw.Layout(width="300px", height="auto")
        # self.raman_tensors_table.layout = ipw.Layout(width="auto", height="auto")

    def _on_site_change(self, change):
        self._model.on_site_selection_change(change["new"])
