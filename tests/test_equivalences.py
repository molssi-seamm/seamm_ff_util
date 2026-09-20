#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for atom types that the equivalence table does not cover.

An equivalence table need not cover every atom type. A combined forcefield takes
its equivalences from one of the forcefields it is built from, while another may
contribute types of its own that have none -- and need none, because it gives an
explicit parameter for every term they appear in.
"""

import pytest

import seamm_ff_util  # noqa: F401


def test_unknown_type_reports_the_missing_term(pcff):
    """A type with no equivalence must not raise a bare KeyError.

    Indexing the equivalence table for such a type raised ``KeyError: '<type>'``,
    which named the type rather than the term and read like the type was the
    problem. The real problem is that the term has no parameters by any route.
    """
    with pytest.raises(RuntimeError, match="No bond parameters for c-nosuchtype"):
        pcff.bond_parameters("c", "nosuchtype")


def test_unknown_type_in_angle(pcff):
    with pytest.raises(RuntimeError) as excinfo:
        pcff.angle_parameters("c", "c", "nosuchtype")

    message = str(excinfo.value)
    assert "angle parameters for c-c-nosuchtype" in message
    assert "No equivalences defined" in message


def test_unknown_type_in_torsion(pcff):
    with pytest.raises(RuntimeError, match="No torsion parameters"):
        pcff.torsion_parameters("c", "c", "c", "nosuchtype")


def test_known_types_still_use_equivalences(pcff):
    """The fallback itself must still work for types the table does cover."""
    assert pcff._have_equivalences("equivalence", "c", "h")
    assert not pcff._have_equivalences("equivalence", "c", "nosuchtype")
    # None is for the terms that take an optional second type.
    assert pcff._have_equivalences("equivalence", "c", None)


def test_assign_forcefield_returns_its_warnings(pcff, configuration):
    """Warnings must come back to the caller, not only go to the log.

    A step needs to put them in its output, where the user will see them. Sent to
    the logger they end up in whatever the job's stdout happens to be -- a SLURM
    output file, or nothing at all.
    """
    configuration.from_smiles("CCO")

    warnings = pcff.assign_forcefield(configuration)

    assert isinstance(warnings, list)
    for warning in warnings:
        assert isinstance(warning, str)


def test_charge_adjustment_is_reported(pcff, configuration):
    """When the charges have to be adjusted, say so, and by how much.

    Ethanol's forcefield charges sum to zero, so telling the configuration it is a
    cation forces the adjustment.
    """
    configuration.from_smiles("CCO")
    configuration.charge = 1

    warnings = pcff.assign_forcefield(configuration)

    assert len(warnings) == 1
    assert "does not match the charge of the system, 1" in warnings[0]
    # The adjustment is one electron over the 9 atoms, and must be reported with
    # enough figures to be meaningful -- it used to be rounded to three decimals,
    # which showed a real adjustment as 0.000.
    assert "0.111111" in warnings[0]
    key = f"charges_{pcff.current_forcefield}"
    assert sum(configuration.atoms[key]) == pytest.approx(1.0)
