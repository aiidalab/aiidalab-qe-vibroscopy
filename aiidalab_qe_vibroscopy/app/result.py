"""Bands results view widgets

"""
from __future__ import annotations


from aiidalab_qe.common.panel import ResultPanel
from aiidalab_qe.common.bandpdoswidget import BandPdosPlotly

import numpy as np

from ..utils.raman.result import export_iramanworkchain_data

from ..utils.phonons.result import export_phononworkchain_data

from ..utils.euphonic import (
    export_euphonic_data,
    EuphonicSuperWidget,
    DowloadYamlHdf5Widget,
)
import plotly.graph_objects as go
import ipywidgets as ipw

from ..utils.raman.result import SpectrumPlotWidget, ActiveModesWidget


class PhononBandPdosPlotly(BandPdosPlotly):
    def __init__(self, bands_data=None, pdos_data=None):
        super().__init__(bands_data, pdos_data)
        self._bands_yaxis = go.layout.YAxis(
            title=dict(text="Phonon Bands (THz)", standoff=1),
            side="left",
            showgrid=True,
            showline=True,
            zeroline=True,
            range=self.SETTINGS["vertical_range_bands"],
            fixedrange=False,
            automargin=True,
            ticks="inside",
            linewidth=2,
            linecolor=self.SETTINGS["axis_linecolor"],
            tickwidth=2,
            zerolinewidth=2,
        )

        paths = self.bands_data.get("paths")
        slider_bands = go.layout.xaxis.Rangeslider(
            thickness=0.08,
            range=[0, paths[-1]["x"][-1]],
        )
        self._bands_xaxis = go.layout.XAxis(
            title="q-points",
            range=[0, paths[-1]["x"][-1]],
            showgrid=True,
            showline=True,
            tickmode="array",
            rangeslider=slider_bands,
            fixedrange=False,
            tickvals=self.bands_data["pathlabels"][1],  # ,self.band_labels[1],
            ticktext=self.bands_data["pathlabels"][0],  # self.band_labels[0],
            showticklabels=True,
            linecolor=self.SETTINGS["axis_linecolor"],
            mirror=True,
            linewidth=2,
            type="linear",
        )


class Result(ResultPanel):
    """
    The idea is that this Panel should be divided in sub panels,
    one for each section of properties: phonons, spectroscopies, inelastic neutron scattering.
    """

    title = "Vibrational Structure"
    workchain_label = "iraman"
    children_result_widget = ()

    def _update_view(self):

        children_result_widget = ()
        tab_titles = []  # this is needed to name the sub panels

        spectra_data = export_iramanworkchain_data(self.node)
        phonon_data = export_phononworkchain_data(self.node)
        ins_data = export_euphonic_data(self.node)

        if phonon_data:

            phonon_children = ()
            if phonon_data["bands"] or phonon_data["pdos"]:
                _bands_plot_view_class = PhononBandPdosPlotly(
                    bands_data=phonon_data["bands"][0],
                    pdos_data=phonon_data["pdos"][0],
                )

                download_widget = DowloadYamlHdf5Widget(
                    phonopy_node=self.node.outputs.vibronic.phonon_pdos.creator
                )

                phonon_children += (
                    _bands_plot_view_class.bandspdosfigure,
                    download_widget,
                )

            if phonon_data["thermo"]:
                import plotly.graph_objects as go

                T = phonon_data["thermo"][0][0]
                F = phonon_data["thermo"][0][1]
                F_units = phonon_data["thermo"][0][2]
                E = phonon_data["thermo"][0][3]
                E_units = phonon_data["thermo"][0][4]
                Cv = phonon_data["thermo"][0][5]
                Cv_units = phonon_data["thermo"][0][6]

                g = go.FigureWidget(
                    layout=go.Layout(
                        title=dict(text="Thermal properties"),
                        barmode="overlay",
                    )
                )
                g.update_layout(
                    xaxis=dict(
                        title="Temperature (K)",
                        linecolor="black",
                        linewidth=2,
                        showline=True,
                    ),
                    yaxis=dict(linecolor="black", linewidth=2, showline=True),
                    plot_bgcolor="white",
                )
                g.add_scatter(x=T, y=F, name=f"Helmoltz Free Energy ({F_units})")
                g.add_scatter(x=T, y=E, name=f"Entropy ({E_units})")
                g.add_scatter(x=T, y=Cv, name=f"Specific Heat-V=const ({Cv_units})")

                phonon_children += (g,)

            tab_titles.append("Phonon properties")

            children_result_widget += (
                ipw.VBox(
                    children=phonon_children,
                    layout=ipw.Layout(
                        width="100%",
                    ),
                ),
            )  # the comma is required! otherwise the tuple is not detected.
        # euphonic
        if ins_data:
            intensity_maps = EuphonicSuperWidget(fc=ins_data["fc"])
            children_result_widget += (intensity_maps,)
            tab_titles.append(f"Inelastic Neutrons")

        if spectra_data:

            # Here we should provide the possibility to have both IR and Raman,
            # as the new logic can provide both at the same time.
            # We are gonna use the same widget, providing the correct spectrum_type: "Raman" or "Ir".
            children_spectra = ()
            for spectrum, data in spectra_data.items():

                if not data:
                    continue

                elif isinstance(data, str):
                    # No Modes are detected. So we explain why
                    no_mode_widget = ipw.HTML(data)
                    explanation_widget = ipw.HTML(
                        "This may be due to the fact that the current implementation of aiida-vibroscopy plugin only considers first-order effects."
                    )

                    children_spectra += (
                        ipw.VBox([no_mode_widget, explanation_widget]),
                    )

                else:
                    subwidget_title = ipw.HTML(f"<h3>{spectrum} spectroscopy</h3>")
                    spectrum_widget = SpectrumPlotWidget(
                        node=self.node, output_node=data, spectrum_type=spectrum
                    )
                    modes_animation = ActiveModesWidget(
                        node=self.node, output_node=data, spectrum_type=spectrum
                    )

                    children_spectra += (
                        ipw.VBox([subwidget_title, spectrum_widget, modes_animation]),
                    )
            children_result_widget += (
                ipw.VBox(
                    children=children_spectra,
                    layout=ipw.Layout(
                        width="100%",
                    ),
                ),
            )
            tab_titles.append(f"Raman/IR spectra")

        self.result_tabs = ipw.Tab(children=children_result_widget)

        for title_index in range(len(tab_titles)):
            self.result_tabs.set_title(title_index, tab_titles[title_index])

        self.children = [self.result_tabs]
