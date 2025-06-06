#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `seamm_ff_util` package."""

import json
import seamm_ff_util  # noqa: F401


def test_angle_angle_torsion_1_explicit(pcff):
    """Test of angle_angle_torsion_1 parameters, which should find
    explicit ones"""

    expected = {
        "Theta0_L": "110.7700",
        "Theta0_R": "108.4000",
        "original Theta0_L": "110.7700",
        "original Theta0_R": "108.4000",
        "reference": "1",
        "version": "1.0",
        "K": "-14.3155",
        "original K": "-14.3155",
    }

    i = "h"
    j = "c"
    k = "c"
    l = "c_0"  # noqa: E741
    ptype, key, form, parameters = pcff.angle_angle_torsion_1_parameters(i, j, k, l)
    assert ptype == "explicit"
    assert key == ("c_0", "c", "c", "h")
    if parameters != expected:
        print(json.dumps(parameters, indent=4))
    assert parameters == expected


def test_angle_angle_torsion_1_explicit_lkji(pcff):
    """known angle_angle_torsion_1 parameters, ordered backwards"""

    expected = {
        "Theta0_L": "108.4000",
        "Theta0_R": "110.7700",
        "original Theta0_L": "108.4000",
        "original Theta0_R": "110.7700",
        "reference": "1",
        "version": "1.0",
        "K": "-14.3155",
        "original K": "-14.3155",
    }

    i = "c_0"
    j = "c"
    k = "c"
    l = "h"  # noqa: E741
    ptype, key, form, parameters = pcff.angle_angle_torsion_1_parameters(i, j, k, l)
    assert ptype == "explicit"
    assert key == ("c_0", "c", "c", "h")
    if parameters != expected:
        print(json.dumps(parameters, indent=4))
    assert parameters == expected


def test_angle_angle_torsion_1_equivalent(pcff):
    """Simple test of angle_angle_torsion_1 parameters using equivalencies"""

    expected = {
        "Theta0_L": "111.0000",
        "Theta0_R": "120.0500",
        "original Theta0_L": "111.0000",
        "original Theta0_R": "120.0500",
        "reference": "1",
        "version": "1.0",
        "K": "-5.8888",
        "original K": "-5.8888",
    }

    i = "h"
    j = "c"
    k = "c5"
    l = "c5"  # noqa: E741
    ptype, key, form, parameters = pcff.angle_angle_torsion_1_parameters(i, j, k, l)
    assert ptype == "equivalent"
    assert key == ("h", "c", "cp", "cp")
    if parameters != expected:
        print(json.dumps(parameters, indent=4))
    assert parameters == expected
