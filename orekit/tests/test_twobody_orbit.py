""" Test cases for two-body orbit behaviors """

# Import required orekit data in current directory
import pytest
from astropy import units as u
import numpy as np
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import assert_allclose
from orekit_orbit import OrekitOrbit
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from utils import setup_orekit_env

# Start orekit virtual machine and load data
vm = setup_orekit_env()


def test_orekit_orbit_from_rv():

    # Initial RV elements
    r = [-6045, -3490, 2500] * u.km
    v = [-3.457, 6.618, 2.533] * u.km / u.s

    # Generate orekit and poliastro orbits
    ss_orekit = OrekitOrbit.from_vectors(Earth, r, v)
    ss_poliastro = Orbit.from_vectors(Earth, r, v)

    # Assert position and velocity vectors
    assert_quantity_allclose(ss_poliastro.a, ss_orekit.a)
    assert_quantity_allclose(ss_poliastro.ecc, ss_orekit.ecc)
    assert_quantity_allclose(ss_poliastro.inc, ss_orekit.inc)
    assert_quantity_allclose(ss_poliastro.raan, ss_orekit.raan)
    assert_quantity_allclose(ss_poliastro.argp, ss_orekit.argp)
    assert_quantity_allclose(ss_poliastro.nu, ss_orekit.nu)


def test_orekit_orbit_from_classical():

    # Initial COE elements
    a = 1.523679 * u.AU
    ecc = 0.093315 * u.one
    inc = 1.85 * u.deg
    raan = 49.562 * u.deg
    argp = 286.537 * u.deg
    nu = 23.33 * u.deg

    # Generate COE element set
    coe = [a, ecc, inc, raan, argp, nu]

    # Generate orekit and poliastro orbits
    ss_orekit = OrekitOrbit.from_classical(Earth, *coe)
    ss_poliastro = Orbit.from_classical(Earth, *coe)

    # Assert position and velocity vectors
    assert_quantity_allclose(ss_poliastro.r, ss_orekit.r)
    assert_quantity_allclose(ss_poliastro.v, ss_orekit.v)
