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
