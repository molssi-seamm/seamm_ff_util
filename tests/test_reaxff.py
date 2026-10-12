#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""ReaxFF: converting a Reax ffield to the MolSSI format keeps every digit."""

from pathlib import Path

from seamm_ff_util.reaxff import (
    _ffield_value,
    read_reaxff_parameters,
    reaxff_to_molssi_format,
)

DATA = Path(__file__).parent / "data"


def _values(text, section, label):
    """The values on the line of ``label`` in ``section`` of a .frc text."""
    lines = text.split(f"#{section}")[1].splitlines()
    for line in lines:
        if line.startswith("#"):
            break
        if label in line.split():
            return line.split()
    raise KeyError(label)


def test_frc_keeps_every_digit():
    """tabulate's default kept 6 significant figures: 153.3934 -> 153.393."""
    parameters = read_reaxff_parameters(DATA / "ffield.reax.h")
    text = reaxff_to_molssi_format(parameters, "test", "2026.10.12")
    assert "9.54691234" in _values(text, "reaxff_general_parameters", "Pboc,2")
    bond = [s for s in text.split("\n#") if s.startswith("reaxff_bond_parameters")]
    assert any("153.3934" in section.split() for section in bond)


def test_ffield_values():
    """LAMMPS's file: four decimals when exact, else every digit."""
    assert _ffield_value("142.2858", 8) == "142.2858"
    assert _ffield_value(-70.1292, 8) == "-70.1292"
    assert _ffield_value(50, 9) == "  50.0000"
    assert _ffield_value(0.12345678, 8) == "0.12345678"
    assert float(_ffield_value(1e-6, 8)) == 1e-6
