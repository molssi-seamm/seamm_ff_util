#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests of rigid bonds and angles, which override the flexible terms."""

import pytest

from seamm_ff_util import Forcefield

FRC = """\
!BIOSYM forcefield 1

#define water

! Version  Ref    Function                         Label
!--------- ---    ------------------------------   ----------
2026.10.15   1    atom_types                       water
2026.10.15   1    equivalence                      water
2026.10.15   1    nonbond(12-6)                    water
2026.10.15   1    quadratic_bond                   water
2026.10.15   1    quadratic_angle                  water

#define water-rigid

! Version  Ref    Function                         Label
!--------- ---    ------------------------------   ----------
2026.10.15   1    atom_types                       water
2026.10.15   1    equivalence                      water
2026.10.15   1    nonbond(12-6)                    water
2026.10.15   1    quadratic_bond                   water
2026.10.15   1    quadratic_angle                  water
2026.10.15   1    rigid_bond                       rigid
2026.10.15   1    rigid_angle                      rigid

#atom_types water

! Version  Ref  Type     Mass      Element     Bonds      Comment
!--------- ---  -----  ----------  -------   ----------   -----------
2026.10.15   1  ow     15.9994      O           2        water O
2026.10.15   1  hw      1.0080      H           1        water H
2026.10.15   1  c       12.011      C           4        carbon

#equivalence water

! Version  Ref  Type     NonB     Bond     Angle    Torsion  OOP
!--------- ---  -------  -------  -------  -------  -------  -------
2026.10.15   1  ow       ow       ow       ow       ow       ow
2026.10.15   1  hw       hw       hw       hw       hw       hw
2026.10.15   1  c        c        c        c        c        c

#nonbond(12-6) water

@type sigma-eps
@combination geometric

! Version  Ref    I        sigma      eps
!--------- ---  -------  ---------  --------
2026.10.15   1  ow       3.166      0.1553
2026.10.15   1  hw       0.0        0.0
2026.10.15   1  c        3.5        0.066

#quadratic_bond water

! Version  Ref     I       J       R0         K2
!--------- ---  ------- -------  -------    --------
2026.10.15   1   hw      ow       1.012      529.6
2026.10.15   1   c       hw       1.09       340.0

#quadratic_angle water

! Version  Ref     I       J       K       Theta0     K2
!--------- ---  ------- ------- -------  -------   --------
2026.10.15   1   hw      ow      hw       113.24     37.95

#rigid_bond rigid

! Version  Ref     I       J       R0
!--------- ---  ------- -------  -------
2026.10.15   1   hw      ow       1.0

#rigid_angle rigid

! Version  Ref     I       J       K       Theta0
!--------- ---  ------- ------- -------  -------
2026.10.15   1   hw      ow      hw       109.47

#end
"""


@pytest.fixture(scope="module")
def frc(tmp_path_factory):
    path = tmp_path_factory.mktemp("ff") / "water.frc"
    path.write_text(FRC)
    return str(path)


def test_flexible(frc):
    ff = Forcefield(frc)
    ff.initialize_biosym_forcefield("water")
    ptype, key, form, values = ff.bond_parameters("ow", "hw")
    assert form == "quadratic_bond"
    assert float(values["R0"]) == pytest.approx(1.012)
    ptype, key, form, values = ff.angle_parameters("hw", "ow", "hw")
    assert form == "quadratic_angle"


def test_rigid_overrides_flexible(frc):
    ff = Forcefield(frc)
    ff.initialize_biosym_forcefield("water-rigid")
    assert ff.ff["terms"]["bond"][0] == "rigid_bond"
    ptype, key, form, values = ff.bond_parameters("ow", "hw")
    assert form == "rigid_bond"
    assert float(values["R0"]) == pytest.approx(1.0)
    ptype, key, form, values = ff.angle_parameters("hw", "ow", "hw")
    assert form == "rigid_angle"
    assert float(values["Theta0"]) == pytest.approx(109.47)


def test_rigid_only_where_given(frc):
    """Other bonds keep their flexible parameters."""
    ff = Forcefield(frc)
    ff.initialize_biosym_forcefield("water-rigid")
    ptype, key, form, values = ff.bond_parameters("c", "hw")
    assert form == "quadratic_bond"
