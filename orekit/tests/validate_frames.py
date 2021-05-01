""" Validate poliastro frames against Orekit ones """

from itertools import product

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import CartesianRepresentation
from astropy.tests.helper import assert_quantity_allclose
from numpy.linalg import norm
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.bodies import CelestialBodyFactory
from org.orekit.frames import Frame, FramesFactory, Transform
from org.orekit.time import AbsoluteDate
from org.orekit.utils import (
    IERSConventions,
    PVCoordinatesProvider,
    TimeStampedPVCoordinates,
)
from poliastro.bodies import Earth as Earth_poliastro
from poliastro.bodies import Jupiter as Jupiter_poliastro
from poliastro.bodies import Mars as Mars_poliastro
from poliastro.bodies import Mercury as Mercury_poliastro
from poliastro.bodies import Moon as Moon_poliastro
from poliastro.bodies import Neptune as Neptune_poliastro
from poliastro.bodies import Saturn as Saturn_poliastro
from poliastro.bodies import Sun as Sun_poliastro
from poliastro.bodies import Uranus as Uranus_poliastro
from poliastro.bodies import Venus as Venus_poliastro
from poliastro.constants import J2000 as J2000_POLIASTRO
from poliastro.frames.equatorial import GCRS as GCRS_poliastro
from poliastro.frames.equatorial import HCRS as HCRS_poliastro
from poliastro.frames.equatorial import JupiterICRS as JupiterICRS_poliastro
from poliastro.frames.equatorial import MarsICRS as MarsICRS_poliastro
from poliastro.frames.equatorial import MercuryICRS as MercuryICRS_poliastro
from poliastro.frames.equatorial import MoonICRS as MoonICRS_poliastro
from poliastro.frames.equatorial import NeptuneICRS as NeptuneICRS_poliastro
from poliastro.frames.equatorial import SaturnICRS as SaturnICRS_poliastro
from poliastro.frames.equatorial import UranusICRS as UranusICRS_poliastro
from poliastro.frames.equatorial import VenusICRS as VenusICRS_poliastro
from poliastro.frames.fixed import ITRS as ITRS_poliastro
from poliastro.frames.fixed import JupiterFixed as JupiterFixed_poliastro
from poliastro.frames.fixed import MarsFixed as MarsFixed_poliastro
from poliastro.frames.fixed import MercuryFixed as MercuryFixed_poliastro
from poliastro.frames.fixed import MoonFixed as MoonFixed_poliastro
from poliastro.frames.fixed import NeptuneFixed as NeptuneFixed_poliastro
from poliastro.frames.fixed import SaturnFixed as SaturnFixed_poliastro
from poliastro.frames.fixed import SunFixed as SunFixed_poliastro
from poliastro.frames.fixed import UranusFixed as UranusFixed_poliastro
from poliastro.frames.fixed import VenusFixed as VenusFixed_poliastro

import orekit
from orekit.pyhelpers import setup_orekit_curdir

# Setup orekit virtual machine and associated data
VM = orekit.initVM()
setup_orekit_curdir("orekit-data.zip")

# All interesting 3D directions
R_SET, V_SET = [list(product([-1, 0, 1], repeat=3)) for _ in range(2)]

# Retrieve celestial bodies from orekit
Sun_orekit = CelestialBodyFactory.getSun()
Mercury_orekit = CelestialBodyFactory.getMercury()
Venus_orekit = CelestialBodyFactory.getVenus()
Earth_orekit = CelestialBodyFactory.getEarth()
Moon_orekit = CelestialBodyFactory.getMoon()
Mars_orekit = CelestialBodyFactory.getMars()
Jupiter_orekit = CelestialBodyFactory.getJupiter()
Saturn_orekit = CelestialBodyFactory.getSaturn()
Uranus_orekit = CelestialBodyFactory.getUranus()
Neptune_orekit = CelestialBodyFactory.getNeptune()

# Name of the bodies
BODIES_NAMES = [
    "Sun",
    "Mercury",
    "Venus",
    "Earth",
    "Moon",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
]

# orekit: bodies, inertial and fixed frames
OREKIT_BODIES = [
    Sun_orekit,
    Mercury_orekit,
    Venus_orekit,
    Earth_orekit,
    Moon_orekit,
    Mars_orekit,
    Jupiter_orekit,
    Saturn_orekit,
    Uranus_orekit,
    Neptune_orekit,
]
OREKIT_FIXED_FRAMES = [body.getBodyOrientedFrame() for body in OREKIT_BODIES]

# poliastro: bodies, inertial and fixed frames
POLIASTRO_BODIES = [
    Sun_poliastro,
    Mercury_poliastro,
    Venus_poliastro,
    Earth_poliastro,
    Moon_poliastro,
    Mars_poliastro,
    Jupiter_poliastro,
    Saturn_poliastro,
    Uranus_poliastro,
    Neptune_poliastro,
]
POLIASTRO_ICRS_FRAMES = [
    HCRS_poliastro,
    MercuryICRS_poliastro,
    VenusICRS_poliastro,
    GCRS_poliastro,
    MoonICRS_poliastro,
    MarsICRS_poliastro,
    JupiterICRS_poliastro,
    SaturnICRS_poliastro,
    UranusICRS_poliastro,
    NeptuneICRS_poliastro,
]
POLIASTRO_FIXED_FRAMES = [
    SunFixed_poliastro,
    MercuryFixed_poliastro,
    VenusFixed_poliastro,
    ITRS_poliastro,
    MoonFixed_poliastro,
    MarsFixed_poliastro,
    JupiterFixed_poliastro,
    SaturnFixed_poliastro,
    UranusFixed_poliastro,
    NeptuneFixed_poliastro,
]


# Collect both API data in two dictionaries
OREKIT_BODIES_AND_FRAMES = dict(
    zip(BODIES_NAMES, zip(OREKIT_BODIES, OREKIT_FIXED_FRAMES))
)
POLIASTRO_BODIES_AND_FRAMES = dict(
    zip(
        BODIES_NAMES,
        zip(POLIASTRO_BODIES, POLIASTRO_ICRS_FRAMES, POLIASTRO_FIXED_FRAMES),
    )
)


# Macros for J2000 and ICRF frame from Orekit API
J2000_OREKIT = AbsoluteDate.J2000_EPOCH
ICRF_FRAME_OREKIT = FramesFactory.getICRF()
GCRF_FRAME_OREKIT = FramesFactory.getGCRF()
ITRF_FRAME_OREKIT = FramesFactory.getITRF(IERSConventions.IERS_2010, False)


# Some of tests are marked as XFAIL since Orekit implements the data from IAU
# WGCCRE 2009 report while poliastro uses IAU WGCCRE 2015 one


@pytest.mark.parametrize("r_vec", R_SET)
@pytest.mark.parametrize("v_vec", V_SET)
@pytest.mark.parametrize(
    "body_name",
    [
        "Sun",
        pytest.param("Mercury", marks=pytest.mark.xfail),
        "Venus",
        "Moon",
        pytest.param("Mars", marks=pytest.mark.xfail),
        "Jupiter",
        "Saturn",
        "Uranus",
        pytest.param("Neptune", marks=pytest.mark.xfail),
    ],
)
def validate_from_body_intertial_to_body_fixed(body_name, r_vec, v_vec):

    # poliastro: collect body information
    (
        BODY_POLIASTRO,
        BODY_ICRF_FRAME_POLIASTRO,
        BODY_FIXED_FRAME_POLIASTRO,
    ) = POLIASTRO_BODIES_AND_FRAMES[body_name]

    # Compute for the norm of position and velocity vectors
    r_norm, v_norm = [norm(vec) for vec in [r_vec, v_vec]]
    R = BODY_POLIASTRO.R.to(u.m).value

    # Correction factor to normalize position and velocity vectors
    k_r = R / r_norm if r_norm != 0 else 1.00
    k_v = 1 / v_norm if v_norm != 0 else 1.00

    # Make a position vector who's norm is equal to the body's radius. Make a
    # unitary velocity vector. Units are in [m] and [m / s].
    rx, ry, rz = [float(k_r * r_i) for r_i in r_vec]
    vx, vy, vz = [float(k_v * v_i) for v_i in v_vec]

    # poliastro: build r_vec and v_vec wrt inertial body frame
    xyz_poliastro = CartesianRepresentation(rx * u.m, ry * u.m, rz * u.m)
    coords_wrt_bodyICRS_poliastro = BODY_ICRF_FRAME_POLIASTRO(xyz_poliastro)

    # poliastro: convert from inertial to fixed frame at given epoch
    coords_wrt_bodyFIXED_poliastro = (
        coords_wrt_bodyICRS_poliastro.transform_to(
            BODY_FIXED_FRAME_POLIASTRO(obstime=J2000_POLIASTRO)
        )
        .represent_as(CartesianRepresentation)
        .xyz
    )

    # orekit: collect body information
    (BODY_OREKIT, BODY_FIXED_FRAME_OREKIT,) = OREKIT_BODIES_AND_FRAMES[body_name]

    # orekit: build r_vec and v_vec wrt inertial body frame
    xyz_orekit = Vector3D(rx, ry, rz)
    uvw_orekit = Vector3D(vx, vy, vz)
    coords_wrt_bodyICRS_orekit = TimeStampedPVCoordinates(
        J2000_OREKIT, xyz_orekit, uvw_orekit
    )

    # orekit: create bodyICRS frame as pure translation of ICRF one
    coords_body_wrt_ICRF_orekit = PVCoordinatesProvider.cast_(
        BODY_OREKIT
    ).getPVCoordinates(J2000_OREKIT, ICRF_FRAME_OREKIT)
    BODY_ICRF_FRAME_OREKIT = Frame(
        ICRF_FRAME_OREKIT,
        Transform(J2000_OREKIT, coords_body_wrt_ICRF_orekit.negate()),
        body_name.capitalize() + "ICRF",
    )

    # orekit: build conversion between BodyICRF and BodyFixed frames
    bodyICRF_to_bodyFIXED_orekit = BODY_ICRF_FRAME_OREKIT.getTransformTo(
        BODY_FIXED_FRAME_OREKIT, J2000_OREKIT,
    )

    # orekit: convert from inertial coordinates to non-inertial ones
    coords_orekit_fixed_raw = (
        bodyICRF_to_bodyFIXED_orekit.transformPVCoordinates(coords_wrt_bodyICRS_orekit)
        .getPosition()
        .toArray()
    )
    coords_wrt_bodyFIXED_orekit = np.asarray(coords_orekit_fixed_raw) * u.m

    # Check position conversion
    assert_quantity_allclose(
        coords_wrt_bodyFIXED_poliastro,
        coords_wrt_bodyFIXED_orekit,
        atol=1e-5 * u.m,
        rtol=1e-7,
    )


@pytest.mark.parametrize("r_vec", R_SET)
@pytest.mark.parametrize("v_vec", V_SET)
def validate_GCRF_to_ITRF(r_vec, v_vec):

    # Compute for the norm of position and velocity vectors
    r_norm, v_norm = [norm(vec) for vec in [r_vec, v_vec]]
    R = Earth_poliastro.R.to(u.m).value

    # Correction factor to normalize position and velocity vectors
    k_r = R / r_norm if r_norm != 0 else 1.00
    k_v = 1 / v_norm if v_norm != 0 else 1.00

    # Make a position vector who's norm is equal to the body's radius. Make a
    # unitary velocity vector. Units are in [m] and [m / s].
    rx, ry, rz = [float(k_r * r_i) for r_i in r_vec]
    vx, vy, vz = [float(k_v * v_i) for v_i in v_vec]

    # orekit: build r_vec and v_vec wrt inertial body frame
    xyz_orekit = Vector3D(rx, ry, rz)
    uvw_orekit = Vector3D(vx, vy, vz)
    coords_GCRF_orekit = TimeStampedPVCoordinates(J2000_OREKIT, xyz_orekit, uvw_orekit)

    # orekit: build conversion between GCRF and ITRF
    GCRF_TO_ITRF_OREKIT = GCRF_FRAME_OREKIT.getTransformTo(
        ITRF_FRAME_OREKIT, J2000_OREKIT,
    )

    # orekit: convert from GCRF to ITRF using previous built conversion
    coords_ITRF_orekit = (
        GCRF_TO_ITRF_OREKIT.transformPVCoordinates(coords_GCRF_orekit)
        .getPosition()
        .toArray()
    )
    coords_ITRF_orekit = np.asarray(coords_ITRF_orekit) * u.m

    # poliastro: build r_vec and v_vec wrt GCRF
    xyz_poliastro = CartesianRepresentation(rx * u.m, ry * u.m, rz * u.m)
    coords_GCRS_poliastro = GCRS_poliastro(xyz_poliastro)

    # poliastro: convert from inertial to fixed frame at given epoch
    coords_ITRS_poliastro = (
        coords_GCRS_poliastro.transform_to(ITRS_poliastro(obstime=J2000_POLIASTRO))
        .represent_as(CartesianRepresentation)
        .xyz
    )

    # Check position conversion
    assert_quantity_allclose(
        coords_ITRS_poliastro, coords_ITRF_orekit, atol=1e-3 * u.m, rtol=1e-2,
    )
