#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Fixtires for testing the 'forcefield' package."""

import pytest
from forcefield.forcefield import Forcefield  # nopep8
from forcefield.ff_assigner import FFAssigner


@pytest.fixture(scope='session')
def pcff():
    """A forcefield object initialized with PCFF
    """
    pcff = Forcefield('data/pcff2017.frc')
    pcff.initialize_biosym_forcefield()
    return pcff


@pytest.fixture(scope='session')
def pcff_assigner(pcff):
    """A forcefield object initialized with PCFF
    """
    pcff_assigner = FFAssigner(pcff)
    return pcff_assigner
