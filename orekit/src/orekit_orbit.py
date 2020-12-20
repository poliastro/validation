""" Simulate Orekit orbit as poliastro one """

import os

import numpy as np
from astropy import units as u
from orekit.pyhelpers import setup_orekit_curdir
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.frames import FramesFactory
from org.orekit.orbits import KeplerianOrbit, OrbitType, PositionAngle
from org.orekit.propagation.analytical import KeplerianPropagator
from org.orekit.time import AbsoluteDate
from org.orekit.utils import Constants, PVCoordinates
from utils import setup_orekit_env

import orekit

# Start orekit virtual machine and load data
vm = setup_orekit_env()


class OrekitOrbit:
    """ Generate poliastro likely orbit from orekit API """

    def __init__(self, attractor, state, epoch, frame):
        """Initializes the OrekitOrbit from a KeplerianOrbit

        Parameters
        ----------
        state: KeplerianOrbit

        """

        # Minimum attributes for describing an orbit in space-time
        self._attractor = attractor
        self._state = state
        self._epoch = epoch
        self._frame = frame

        # Custom keplerian propagator
        self._propagator = KeplerianPropagator(self._state)

    @property
    def attractor(self):
        """ Attractor of the orbit """
        return self._attractor

    @property
    def epoch(self):
        """ Epoch of the orbit """
        return self._epoch

    @property
    def frame(self):
        """ Frame of the orbit """
        return self._frame

    @property
    def r(self):
        """ Position vector """
        return self._state.getPVCoordinates().position.toArray() * u.m

    @property
    def v(self):
        """ velocity vector """
        return self._state.getPVCoordinates().velocity.toArray() * u.m / u.s

    @property
    def a(self):
        """ Semi-major axis """
        return self._state.getA() * u.m

    @property
    def ecc(self):
        """ Eccentricity """
        return self._state.getE() * u.one

    @property
    def inc(self):
        """ Inclination """
        return self._state.getI() * u.rad

    @property
    def raan(self):
        """ Argument of periapsis """
        # Orekit bounds RAAN between [-pi,pi]. If negative, +2pi is required
        raan = self._state.getRightAscensionOfAscendingNode()

        if -np.pi < raan < 0:
            return (raan + 2.0 * np.pi) * u.rad
        else:
            return raan * u.rad

    @property
    def argp(self):
        """ Argument of perigee """
        return self._state.getPerigeeArgument() * u.rad

    @property
    def nu(self):
        """ True anomaly """
        return self._state.getTrueAnomaly() * u.rad

    @classmethod
    def from_vectors(
        cls,
        attractor,
        r,
        v,
        epoch=AbsoluteDate.J2000_EPOCH,
        frame=FramesFactory.getGCRF(),
    ):
        """ Builds an OrekitOrbit from position and velocity vectors """

        # Get the value of the coordinates for position and velocity vectors
        rx, ry, rz = [float(r_coord.to(u.m).value) for r_coord in r]
        vx, vy, vz = [float(v_coord.to(u.m / u.s).value) for v_coord in v]

        # Convert to Orekit three-dimensional vectors and solve PV instance
        r_vec, v_vec = Vector3D(rx, ry, rz), Vector3D(vx, vy, vz)
        rv_state = PVCoordinates(r_vec, v_vec)

        # Get the gravitational parameter
        k = float(attractor.k.to(u.m ** 3 / u.s ** 2).value)

        # Build a KeplerianOrbit from vectors
        state = KeplerianOrbit(rv_state, frame, epoch, k)

        return cls(attractor, state, epoch, frame)

    @classmethod
    def from_classical(
        cls,
        attractor,
        a,
        ecc,
        inc,
        raan,
        argp,
        nu,
        epoch=AbsoluteDate.J2000_EPOCH,
        frame=FramesFactory.getGCRF(),
    ):

        # Get rid of units and force float precision
        k = float(attractor.k.to(u.m ** 3 / u.s ** 2).value)
        a = float(a.to(u.m).value)
        ecc = float(ecc.to(u.one).value)
        inc = float(inc.to(u.rad).value)
        raan = float(raan.to(u.rad).value)
        argp = float(argp.to(u.rad).value)
        nu = float(nu.to(u.rad).value)

        # Generate a Orekit KeplerianOrbit
        state = KeplerianOrbit(
            a,
            ecc,
            inc,
            argp,
            raan,
            nu,
            PositionAngle.TRUE,
            frame,
            epoch,
            Constants.WGS84_EARTH_MU,
        )

        return cls(attractor, state, epoch, frame)

    def classical(self):
        """ Returns orbit classical elements """
        coe = (
            self.a,
            self.ecc,
            self.inc,
            self.raan,
            self.argp,
            self.nu,
        )

        return coe

    def rv(self):
        """ Returns position and velocity vectors """
        return self.r, self.v

    def propagate(self, tof):
        """ Propagates the orbit a given amount of time """

        # Ensure seconds conversion and floating value
        tof = float(tof.to(u.s).value)

        # Generate the new epoch by shifting orbit's original one
        new_epoch = self.epoch.shiftedBy(tof)

        # Propagate the orbit and get the new position and velocity vectors
        rv_new = self._propagator.propagate(self.epoch, new_epoch).getPVCoordinates()
        r, v = (rv_new.getPosition().toArray()) * u.m, (
            rv_new.getVelocity().toArray()
        ) * u.m / u.s

        # Return a new orbit from
        return OrekitOrbit.from_vectors(self.attractor, r, v, new_epoch, self.frame)
