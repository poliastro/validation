""" Validate poliastro's two-body element conversion against Orekit"""

import numpy as np
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.frames import FramesFactory
from org.orekit.orbits import KeplerianOrbit, PositionAngle
from org.orekit.time import AbsoluteDate
from org.orekit.utils import Constants as C
from org.orekit.utils import PVCoordinates
from poliastro.bodies import Earth
from poliastro.twobody import Orbit

import orekit
from orekit.pyhelpers import setup_orekit_curdir

# Setup orekit virtual machine and associated data
VM = orekit.initVM()
setup_orekit_curdir("orekit-data.zip")


def validate_rv2coe():
    # Initial orbit state vectors, final radius and time of flight
    rx_0, ry_0, rz_0 = [float(ri_0) for ri_0 in [-6045, -3490, 2500]]
    vx_0, vy_0, vz_0 = [float(vi_0) for vi_0 in [-3.457, 6.618, 2.533]]

    # Build the initial Orekit orbit
    k = C.IAU_2015_NOMINAL_EARTH_GM
    r0_vec = Vector3D(rx_0, ry_0, rz_0)
    v0_vec = Vector3D(vx_0, vy_0, vz_0)
    rv_0 = PVCoordinates(r0_vec, v0_vec)
    epoch_00 = AbsoluteDate.J2000_EPOCH
    gcrf_frame = FramesFactory.getGCRF()
    ss0_orekit = KeplerianOrbit(rv_0, gcrf_frame, epoch_00, k)

    # Build poliastro orbit
    r0_vec = [rx_0, ry_0, rz_0] * u.m
    v0_vec = [vx_0, vy_0, vz_0] * u.m / u.s
    ss0_poliastro = Orbit.from_vectors(Earth, r0_vec, v0_vec)

    # Orekit bounds COE angular magnitudes between [-pi, +pi]. Let us define a
    # converter function to map values between [0, +2pi]
    def unbound_angle(angle):
        """ Bound angle between [0, +2pi] """
        if -np.pi <= angle <= 0.0:
            angle += 2 * np.pi

        return angle

    # Map angles
    orekit_raan = unbound_angle(ss0_orekit.getRightAscensionOfAscendingNode())
    orekit_argp = unbound_angle(ss0_orekit.getPerigeeArgument())
    orekit_nu = unbound_angle(ss0_orekit.getTrueAnomaly())

    # Assert classical orbital elements
    assert_quantity_allclose(ss0_poliastro.a, ss0_orekit.getA() * u.m, rtol=1e-6)
    assert_quantity_allclose(ss0_poliastro.ecc, ss0_orekit.getE() * u.one, rtol=1e-6)
    assert_quantity_allclose(ss0_poliastro.inc, ss0_orekit.getI() * u.rad, rtol=1e-6)
    assert_quantity_allclose(ss0_poliastro.raan, orekit_raan * u.rad, rtol=1e-6)
    assert_quantity_allclose(ss0_poliastro.argp, orekit_argp * u.rad, rtol=1e-6)
    assert_quantity_allclose(ss0_poliastro.nu, orekit_nu * u.rad, rtol=1e-6)


def validate_coe2rv():

    # Initial COE of the orbit
    a, ecc, inc, raan, argp, nu = [87800e3, 0.1712, 2.6738, 4.4558, 0.3612, 0.4965]

    # Build Orekit orbit
    k = C.IAU_2015_NOMINAL_EARTH_GM
    epoch_00 = AbsoluteDate.J2000_EPOCH
    gcrf_frame = FramesFactory.getGCRF()
    ss0_orekit = KeplerianOrbit(
        a, ecc, inc, argp, raan, nu, PositionAngle.TRUE, gcrf_frame, epoch_00, k
    )

    # Build poliastro orbit
    ss0_poliastro = Orbit.from_classical(
        Earth,
        a * u.m,
        ecc * u.one,
        inc * u.rad,
        raan * u.rad,
        argp * u.rad,
        nu * u.rad,
    )

    # Retrieve Orekit final state vectors
    r_orekit = ss0_orekit.getPVCoordinates().position.toArray() * u.m
    v_orekit = ss0_orekit.getPVCoordinates().velocity.toArray() * u.m / u.s

    # Retrieve poliastro final state vectors
    r_poliastro, v_poliastro = ss0_poliastro.rv()

    # Assert final state vector
    assert_quantity_allclose(r_poliastro, r_orekit)
    assert_quantity_allclose(v_poliastro, v_orekit)
