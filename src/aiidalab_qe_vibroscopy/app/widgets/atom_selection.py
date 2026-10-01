"""Shared atom-selection controls for vibrational result panels."""

import ipywidgets as ipw


class AtomSelectionControls(ipw.VBox):
    def __init__(self, number_of_atoms, optical=False):
        self.selection = ipw.Text(
            description="Atoms:",
            placeholder="1, 3; 6..8",
            continuous_update=False,
            style={"description_width": "initial"},
            layout=ipw.Layout(width="450px"),
        )
        self.selected_only = ipw.Checkbox(
            description="Show selected curve only", indent=False
        )
        text = (
            f"Use atom numbers 1–{number_of_atoms} in input-structure order (1-based). "
            "Spaces, commas, semicolons and inclusive ranges (6..8) are accepted. "
            "Leave empty to show the original spectrum."
        )
        if optical:
            text += (
                "<br><b>Selected-atom mode participation:</b> mass-weighted motion "
                "in each mode, shown on the total spectrum’s scale. This is not "
                "an additive decomposition of atomic IR/Raman intensities. "
                "Degenerate modes share an averaged participation."
            )
        self.help = ipw.HTML(text)
        self.message = ipw.HTML()
        super().__init__(
            [
                self.selection,
                self.selected_only,
                self.help,
                self.message,
            ]
        )

    def close(self):
        for child in getattr(self, "children", ()):
            child.close()
        super().close()
