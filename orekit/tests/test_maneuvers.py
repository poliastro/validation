""" Validates poliastro's maneuvers against Orekit ones """

import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from orekit_maneuvers import OrekitManeuver
from orekit_orbit import OrekitOrbit
from poliastro.bodies import Earth
from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit

# A collection of APIs to be validated against each other
API_set = {"Orekit": [OrekitOrbit, OrekitManeuver], "poliastro": [Orbit, Maneuver]}


@pytest.mark.parametrize("API", API_set)
def test_3D_hohmann_maneuver(API):
    # Initial orbit data
    r = [7200, -1000, 0] * u.km
    v = [0, 8, 0] * u.km / u.s

    # Expected final state vector computed using GMAT-2020a
    tof_expected = 19790.48929209821 * u.s
    r_expected = [
        -18447.18133044379e3,
        -30659.53028128713e3,
        0.002151253568205197e3,
    ] * u.m
    v_expected = (
        [2.859891380100071e3, -1.720738616179917e3, -2.579688114910805e-04] * u.m / u.s
    )

    # Expected impulses as seen from the PQW local orbit frame
    expected_dVa = [1.47424834733, 0, 0] * u.km / u.s
    expected_dVb = [1.44567622075, 0, 0] * u.km / u.s

    # Unpack API utilities
    OrbitAPI, ManeuverAPI = API_set[API]

    # Initial orbit, maneuver and final orbit
    ss_0 = OrbitAPI.from_vectors(Earth, r, v)
    man = ManeuverAPI.hohmann(ss_0, 35781.34857 * u.km)

    # Orekit is triggered by internal events while poliastro requires
    # propagation after maneuver has been applied
    if API == "Orekit":
        ss_f = ss_0.apply_maneuver(man, tof_expected)
        assert_quantity_allclose(man._dvs[0], expected_dVa, rtol=1e-5)
        assert_quantity_allclose(man._dvs[-1], expected_dVb, rtol=1e-5)
    else:
        ss_f = ss_0.apply_maneuver(man).propagate(1 * u.hour)

    assert_quantity_allclose(ss_f.r, r_expected, atol=1 * u.km)
    assert_quantity_allclose(ss_f.v, v_expected, atol=0.01 * u.km / u.s)
