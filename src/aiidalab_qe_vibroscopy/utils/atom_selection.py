"""Parse the one-based atom selections used by the spectrum viewers."""

import re


def parse_atom_selection(text: str, number_of_atoms: int) -> list[int]:
    """Return sorted, unique zero-based indices, rejecting ambiguous input.

    Like the EMPA atom-index parser, accept whitespace, commas, semicolons
    and inclusive '..' ranges, including whitespace around the range.
    Also accept a surrounding pair of parentheses/brackets and a trailing
    full stop after an integer. Decimal numbers are deliberately rejected.
    Bounds are checked before expanding a range.
    """
    value = text.strip()
    if not value:
        return []
    if value[0] in "([":
        closing = {"(": ")", "[": "]"}[value[0]]
        if not value.endswith(closing):
            raise ValueError("Close the atom list with the matching bracket.")
        value = value[1:-1].strip()
    if not value:
        return []
    # A full stop followed by a separator is punctuation; 5.8 is not two atoms.
    value = re.sub(r"(?<=\d)\.(?=\s|[,;]|$)", "", value)
    value = re.sub(r"\s*\.\.\s*", "..", value)
    value = re.sub(r"[,;]+", " ", value)
    selected = set()
    for token in value.split():
        if not re.fullmatch(r"\d+(?:\.\.\d+)?", token):
            raise ValueError(
                f"Invalid atom selection: {token!r}. Use indices such as 1, 3; 6..8."
            )
        limits = [int(part) for part in token.split("..")]
        first, last = limits[0], limits[-1]
        if first > last:
            raise ValueError(f"Atom range {token!r} is reversed.")
        if first < 1 or last > number_of_atoms:
            raise ValueError(f"Atom indices must be between 1 and {number_of_atoms}.")
        selected.update(range(first - 1, last))
    if not selected:
        raise ValueError("Enter atom indices, or leave the field empty.")
    return sorted(selected)
