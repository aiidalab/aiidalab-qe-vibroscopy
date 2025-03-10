from aiidalab_qe_vibroscopy.utils.euphonic.data.phonopy_interface import (
    generate_force_constant_from_phonopy,
)


def export_euphonic_data(output_vibronic, fermi_energy=None):
    if "phonon_bands" not in output_vibronic:
        return None

    output_set = output_vibronic.phonon_bands

    if any(not element for element in output_set.creator.caller.inputs.structure.pbc):
        vibro_bands = output_set.creator.caller.inputs.phonopy_bands_dict.get_dict()
        # Group the band and band_labels
        band = vibro_bands["band"]
        band_labels = vibro_bands["band_labels"]

        grouped_bands = [
            item
            for sublist in [band_labels[i : i + 2] for i in range(len(band_labels) - 1)]
            for item in sublist
        ]
        grouped_q = [
            [tuple(band[i : i + 3]), tuple(band[i + 3 : i + 6])]
            for i in range(0, len(band) - 3, 3)
        ]
        q_path = {
            "coordinates": grouped_q,
            "labels": grouped_bands,
            "delta_q": 0.01,  # 1/A
        }
    else:
        q_path = None

    phonopy_calc = output_set.creator
    fc = generate_force_constant_from_phonopy(
        phonopy_calc,
        use_euphonic_full_parser=True,
    )
    # bands = compute_bands(fc)
    # pdos = compute_pdos(fc)
    return {
        "fc": fc,
        "q_path": q_path,
    }  # "bands": bands, "pdos": pdos, "thermal": None}
