"""Bands results view widgets

"""



from aiidalab_qe.common.panel import ResultPanel

import numpy as np


def export_phononworkchain_data(node, fermi_energy=None):

    """
    We have multiple choices: BANDS, DOS, THERMODYNAMIC, FORCES.
    """

    import json

    from monty.json import jsanitize

    parameters = {}

    if not "vibronic" in node.outputs:
        return None
    else:
        if (
            not "harmonic" in node.outputs.vibronic
            and not "phonon" in node.outputs.vibronic
        ):
            return None

    output_set = (
        node.outputs.vibronic.harmonic
        if "harmonic" in node.outputs.vibronic
        else node.outputs.vibronic.phonon
    )

    if "output_phonopy" in output_set:
        if "phonon_bands" in output_set.output_phonopy:
            data = json.loads(
                output_set.output_phonopy.phonon_bands._exportcontent(
                    "json", comments=False
                )[0]
            )
            # The fermi energy from band calculation is not robust.
            """data["fermi_level"] = (
                fermi_energy or node.outputs.phonons.band_parameters["fermi_energy"]
            )"""
            # to be optimized: use the above results!!!
            bands = output_set.output_phonopy.phonon_bands.get_bands()
            data["fermi_level"] = 0
            data["Y_label"] = "Dispersion (THz)"

            # it does work now.
            parameters["energy_range"] = {
                "ymin": np.min(bands) - 0.1,
                "ymax": np.max(bands) + 0.1,
            }

            # TODO: THERMOD, FORCES; minors: bands-labels, done: no-fermi-in-dos.

            return [jsanitize(data), parameters, "bands"]
        elif "total_phonon_dos" in output_set.output_phonopy:
            (
                what,
                energy_dos,
                units_omega,
            ) = output_set.output_phonopy.total_phonon_dos.get_x()
            (
                dos_name,
                dos_data,
                units_dos,
            ) = output_set.output_phonopy.total_phonon_dos.get_y()[0]
            dos = []
            # The total dos parsed
            tdos = {
                "label": "Total DOS",
                "x": energy_dos.tolist(),
                "y": dos_data.tolist(),
                "borderColor": "#8A8A8A",  # dark gray
                "backgroundColor": "#8A8A8A",  # light gray
                "backgroundAlpha": "40%",
                "lineStyle": "solid",
            }
            dos.append(tdos)

            parameters["energy_range"] = {
                "ymin": np.min(energy_dos) - 0.1,
                "ymax": np.max(energy_dos) + 0.1,
            }

            data_dict = {
                "fermi_energy": 0,  # I do not want it in my plot
                "dos": dos,
            }

            return [json.loads(json.dumps(data_dict)), parameters, "dos"]
        elif "thermal_properties" in output_set.output_phonopy:
            what, T, units_k = output_set.output_phonopy.thermal_properties.get_x()
            (
                F_name,
                F_data,
                units_F,
            ) = output_set.output_phonopy.thermal_properties.get_y()[0]
            (
                Entropy_name,
                Entropy_data,
                units_entropy,
            ) = output_set.output_phonopy.thermal_properties.get_y()[1]
            (
                Cv_name,
                Cv_data,
                units_Cv,
            ) = output_set.output_phonopy.thermal_properties.get_y()[2]

            return (
                [T, F_data, units_F, Entropy_data, units_entropy, Cv_data, units_Cv],
                [],
                "thermal",
            )

    else:
        return None


class Result(ResultPanel):

    title = "Phonon property"
    workchain_label = "phonons"

    def _update_view(self):
        bands_data = export_phononworkchain_data(self.node)

        if bands_data[2] == "bands":
            _bands_plot_view = BandsPlotWidget(
                bands=[bands_data[0]],
                **bands_data[1],
            )
            self.children = [
                _bands_plot_view,
            ]
        elif bands_data[2] == "dos":
            _bands_plot_view = BandsPlotWidget(
                dos=bands_data[0],
                plot_fermilevel=False,
                show_legend=False,
                **bands_data[1],
            )
            self.children = [
                _bands_plot_view,
            ]

        elif bands_data[2] == "thermal":
            import plotly.graph_objects as go

            T = bands_data[0][0]
            F = bands_data[0][1]
            F_units = bands_data[0][2]
            E = bands_data[0][3]
            E_units = bands_data[0][4]
            Cv = bands_data[0][5]
            Cv_units = bands_data[0][6]

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

            self.children = [
                g,
            ]
