"""Check atom selections against the user-facing indexing convention."""

import pytest

from aiidalab_qe_vibroscopy.utils.atom_selection import parse_atom_selection


@pytest.mark.parametrize(
    "text, expected",
    [
        ("", []),
        ("  ", []),
        ("1 3 6..8", [0, 2, 5, 6, 7]),
        ("1,3 ;5. 8", [0, 2, 4, 7]),
        ("( 1,3 ;5. 8 10 .. 30 44)", [0, 2, 4, 7, *range(9, 30), 43]),
        ("[1; 3 .. 5]", [0, 2, 3, 4]),
        ("1 1, 2..4 3", [0, 1, 2, 3]),
        ("74", [73]),
        ("1..74", list(range(74))),
        ("5.", [4]),
    ],
)
def test_atom_selection(text, expected):
    assert parse_atom_selection(text, 74) == expected


@pytest.mark.parametrize(
    "text",
    [
        "0",
        "-1",
        "75",
        "8..6",
        "1..100000000000",
        "5.8",
        "1...4",
        "1..",
        "..4",
        "1 two 3",
        "(1 3]",
        "(1 3",
        "1 3)",
        "1e1",
        "1/2",
        "1-3",
        "1; __import__('os')",
        ",;",
    ],
)
def test_invalid_atom_selection(text):
    with pytest.raises(ValueError):
        parse_atom_selection(text, 74)
