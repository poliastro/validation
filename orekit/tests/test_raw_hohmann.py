""" Test 3D Hohmann maneuvers from Orekit against poliastro ones """

import astropy.units as u
import numpy as np
from astropy.tests.helper import assert_quantity_allclose
from poliastro.bodies import Earth
from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit

# Setup the virtual machine and data location
import orekit 
from orekit.pyhelpers import setup_orekit_curdir

vm = orekit.initVM() 

# Data orekit-data.zip must be located at same dir/ level that this script
setup_orekit_curdir("src/data/orekit-data.zip") 


# Import different utilities from Orekit API
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.attitudes import LofOffset
from org.orekit.forces.maneuvers import ImpulseManeuver
from org.orekit.frames import FramesFactory, LOFType
from org.orekit.orbits import KeplerianOrbit
from org.orekit.propagation.analytical import KeplerianPropagator
from org.orekit.propagation.events import ApsideDetector
from org.orekit.propagation.events.handlers import StopOnDecreasing
from org.orekit.time import AbsoluteDate
from org.orekit.utils import Constants as C
from org.orekit.utils import PVCoordinates


def test_poliastro_3Dhohmann():

    # Initial orbit state vectors, final radius and time of flight
    rx_0, ry_0, rz_0 = [float(ri_0) for ri_0 in [7200e3, -1000e3, 0]]
    vx_0, vy_0, vz_0 = [float(vi_0) for vi_0 in [0, 8e3, 0]]
    rf_norm = float(35781.34857e3)
    tof = float(19790.48929209821)

    # Build the initial Orekit orbit
    k = C.IAU_2015_NOMINAL_EARTH_GM
    r0_vec = Vector3D(rx_0, ry_0, rz_0)
    v0_vec = Vector3D(vx_0, vy_0, vz_0)
    rv_0 = PVCoordinates(r0_vec, v0_vec)
    epoch_00 = AbsoluteDate.J2000_EPOCH
    gcrf_frame = FramesFactory.getGCRF()
    ss0_orekit = KeplerianOrbit(rv_0, gcrf_frame, epoch_00, k)

    # Final desired orbit radius and auxiliary variables
    r0_norm = ss0_orekit.getA() * (1 - ss0_orekit.getE())
    v0_norm = np.sqrt(2 * k / r0_norm - k / ss0_orekit.getA())
    a_trans = (r0_norm + rf_norm) / 2

    # Compute the magnitude of Hohmann's deltaV
    dv_a = np.sqrt(2 * k / r0_norm - k / a_trans) - v0_norm
    dv_b = np.sqrt(k / rf_norm) - np.sqrt(2 * k / rf_norm - k / a_trans)
    dVa_vec, dVb_vec = [Vector3D(float(dV), float(0), float(0)) for dV in (dv_a, dv_b)]

    # Local orbit frame: X-axis aligned with velocity, Z-axis towards momentum
    attitude_provider = LofOffset(gcrf_frame, LOFType.TNW)

    # Setup the triggers for the impulses execution
    at_periapsis = ApsideDetector(ss0_orekit)
    at_apoapsis = ApsideDetector(ss0_orekit).withHandler(StopOnDecreasing())

    # Build the impulsive maneuvers
    ISP = float(300)
    imp_A = ImpulseManeuver(at_periapsis, attitude_provider, dVa_vec, ISP)
    imp_B = ImpulseManeuver(at_apoapsis, attitude_provider, dVb_vec, ISP)

    # Generate the propagator and add the maneuvers
    propagator = KeplerianPropagator(ss0_orekit, attitude_provider)
    propagator.addEventDetector(imp_A)
    propagator.addEventDetector(imp_B)

    # Expected time of flight
    epoch_ff = epoch_00.shiftedBy(tof)

    # Propagate the Orekit orbit
    rv_f = propagator.propagate(epoch_ff).getPVCoordinates(gcrf_frame)
    rf_orekit, vf_orekit = (
        (list(rv_f.getPosition().toArray()) * u.m),
        (list(rv_f.getVelocity().toArray()) * u.m / u.s),
    )

    # Build the orbit with poliastro
    ss0_poliastro = Orbit.from_vectors(
        Earth, r0_vec.toArray() * u.m, v0_vec.toArray() * u.m / u.s
    )
    man_poliastro = Maneuver.hohmann(ss0_poliastro, rf_norm * u.m)
    rf_poliastro, vf_poliastro = (
        ss0_poliastro.apply_maneuver(man_poliastro).propagate(1 * u.h).rv()
    )

    # Check if final vectors match expected ones predicted by GMAT-NASA 2020a
    rf_expected = [
        -18447.18133044379e3,
        -30659.53028128713e3,
        0.002151253568205197e3,
    ] * u.m
    vf_expected = (
        [2.859891380100071e3, -1.720738616179917e3, -2.579688114910805e-04] * u.m / u.s
    )

    # Test Orekit against GMAT
    assert_quantity_allclose(rf_orekit, rf_expected, atol=50 * u.m, rtol=1e-4)
    assert_quantity_allclose(vf_orekit, vf_expected, atol=0.0035 * u.m / u.s, rtol=1e-5)

    # test poliastro against Orekit
    assert_quantity_allclose(rf_poliastro, rf_orekit, rtol=1e-6)
    assert_quantity_allclose(vf_poliastro, vf_orekit, rtol=1e-6)

    # Assert total time of flight = time of maneuver + 1 hour of propagation
    assert_quantity_allclose(man_poliastro.get_total_time() + 1 * u.h, tof * u.s)
