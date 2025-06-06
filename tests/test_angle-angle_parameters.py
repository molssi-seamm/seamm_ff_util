#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `seamm_ff_util` package."""

import json
import seamm_ff_util  # noqa: F401


def test_angle_angle_explicit(pcff):
    """Test of bond_bond parameters, which should fine explicit ones"""

    expected = {
        "Theta10": "107.6600",
        "Theta20": "105.8500",
        "original Theta10": "107.6600",
        "original Theta20": "105.8500",
        "reference": "1",
        "version": "1.0",
        "K": "3.5475",
        "original K": "3.5475",
    }

    i = "h"
    j = "c"
    k = "h"
    l = "n+"  # noqa: E741
    ptype, key, form, parameters = pcff.angle_angle_parameters(i, j, k, l)
    assert ptype == "explicit"
    assert key == ("h", "c", "h", "n+")
    if parameters != expected:
        print(json.dumps(parameters, indent=4))
    assert parameters == expected


def test_angle_angle_explicit_ljki(pcff):
    """known bond_bond parameters, ordered backwards"""
    i = "h"
    j = "c"
    k = "h"
    l = "n+"  # noqa: E741
    ptype, key, form, parameters = pcff.angle_angle_parameters(i, j, k, l)
    ptype2, key2, form, parameters2 = pcff.angle_angle_parameters(l, j, k, i)
    assert ptype2 == "explicit"
    assert key2 == ("n+", "c", "h", "h")
    assert parameters == parameters2


def test_angle_angle_equivalent(pcff):
    """Simple test of bond_bond parameters using equivalencies"""
    expected = {
        "Theta10": "116.0640",
        "Theta20": "116.0640",
        "original Theta10": "116.0640",
        "original Theta20": "116.0640",
        "reference": "6",
        "version": "2.1",
        "K": "5.9863",
        "original K": "5.9863",
    }

    i = "c5"
    j = "cp"
    k = "c_1"
    l = "c5"  # noqa: E741
    ptype, key, form, parameters = pcff.angle_angle_parameters(i, j, k, l)
    assert ptype == "equivalent"
    assert key == ("cp", "cp", "c_1", "cp")
    if parameters != expected:
        print(json.dumps(parameters, indent=4))
    assert parameters == expected
