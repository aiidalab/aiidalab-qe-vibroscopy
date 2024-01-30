"""Bands results view widgets

"""
from __future__ import annotations


from widget_bandsplot import BandsPlotWidget

from aiidalab_qe.common.panel import ResultPanel

import numpy as np

from ..utils.raman.result import export_iramanworkchain_data

# from ..utils.harmonic.result import export_phononworkchain_data

from ..utils.euphonic.euphonic_widgets import (
    export_phononworkchain_data,
    IntensityFullWidget,
    PowderFullWidget,
)

import ipywidgets as ipw
from ..utils.raman.result import SpectrumPlotWidget, ActiveModesWidget


class Result(ResultPanel):

    title = "Vibrational Structure"
    workchain_label = "iraman"
    children_result_widget = ()

    def _update_view(self):

        children_result_widget = ()

        spectra_data = export_iramanworkchain_data(self.node)
        phonon_data = export_phononworkchain_data(self.node)

        if phonon_data:

            # euphonic
            if phonon_data["bands"]:
                _bands_plot_view = BandsPlotWidget(
                    bands=[phonon_data["bands"][0]],
                    **phonon_data["bands"][1],
                )
                children_result_widget += (_bands_plot_view,)

            # euphonic
            if phonon_data["pdos"]:
                _pdos_plot_view = BandsPlotWidget(
                    dos=phonon_data["pdos"][0],
                    plot_fermilevel=False,
                    show_legend=False,
                    **phonon_data["pdos"][1],
                )
                children_result_widget += (_pdos_plot_view,)

            if phonon_data["thermal"]:
                import plotly.graph_objects as go

                T = phonon_data[0][0]
                F = phonon_data[0][1]
                F_units = phonon_data[0][2]
                E = phonon_data[0][3]
                E_units = phonon_data[0][4]
                Cv = phonon_data[0][5]
                Cv_units = phonon_data[0][6]

                g = go.FigureWidget(
                    layout=go.Layout(
                        title=dict(text="Thermal properties"),
                        barmode="overlay",
                    )
                )
                g.layout.xaxis.title = "Temperature (K)"
                g.add_scatter(x=T, y=F, name=f"Helmoltz Free Energy ({F_units})")
                g.add_scatter(x=T, y=E, name=f"Entropy ({E_units})")
                g.add_scatter(x=T, y=Cv, name=f"Specific Heat-V=const ({Cv_units})")

                children_result_widget += (g,)

            # euphonic
            if phonon_data["fc"]:
                intensity_map = IntensityFullWidget(fc=phonon_data["fc"])
                powder_map = PowderFullWidget(fc=phonon_data["fc"])
                children_result_widget += (intensity_map, powder_map)

        if spectra_data:

            # Here we should provide the possibility to have both IR and Raman,
            # as the new logic can provide both at the same time.
            # We are gonna use the same widget, providing the correct spectrum_type: "Raman" or "Ir".

            for spectrum, data in spectra_data.items():

                if not data:
                    continue

                elif isinstance(data, str):
                    # No Modes are detected. So we explain why
                    no_mode_widget = ipw.HTML(data)
                    explanation_widget = ipw.HTML(
                        "This may be due to the fact that the current implementation of aiida-vibroscopy plugin only considers first-order effects."
                    )

                    children_result_widget += (
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

                    children_result_widget += (
                        ipw.VBox([subwidget_title, spectrum_widget, modes_animation]),
                    )

        self.children = children_result_widget
