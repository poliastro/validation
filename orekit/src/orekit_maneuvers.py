""" Simulate Orekit maneuvers as poliastro ones """

import numpy as np
from astropy import units as u

from utils import vm  # noqa: F401 # isort: split
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.attitudes import LofOffset
from org.orekit.forces.maneuvers import ImpulseManeuver
from org.orekit.frames import LOFType
from org.orekit.propagation.events import ApsideDetector
from org.orekit.propagation.events.handlers import StopOnDecreasing


class OrekitManeuver:
    """ Class for modeling impulsive maneuvers """

    def __init__(self, burns, attitude_provider):
        """ Initializes the maneuver """

        # Save the ImpulsiveManeuver instances and the attitude
        self.burns = burns
        self.attitude_provider = attitude_provider

        # Allocate the final dV vectors
        self._dvs = [dV.getDeltaVSat().toArray() * u.m / u.s for dV in self.burns]

    @classmethod
    def hohmann(cls, orekit_ss, rf_norm, ISP=float(300)):
        """ Computes required impulses for a Hohmann transfer """

        # Get required parameters and auxiliary ones
        k = orekit_ss.attractor.k
        r0_norm = orekit_ss.r_p
        v0_norm = np.sqrt(2 * k / r0_norm - k / orekit_ss.a)
        a_trans = (r0_norm + rf_norm) / 2

        # Compute the magnitude of Hohmann's deltaV
        dv_a = np.sqrt(2 * k / r0_norm - k / a_trans) - v0_norm
        dv_b = np.sqrt(k / rf_norm) - np.sqrt(2 * k / rf_norm - k / a_trans)

        # Local orbit frame: X-axis aligned with velocity, Z-axis towards momentum
        attitude_provider = LofOffset(orekit_ss.frame, LOFType.TNW)

        # Convert previous magnitudes into vectors within the TNW frame, which
        # is aligned with local velocity vector
        dVa_vec, dVb_vec = [
            Vector3D(float(dV.to(u.m / u.s).value), float(0), float(0))
            for dV in (dv_a, dv_b)
        ]

        # Define the events when that trigger the impulses
        at_periapsis = ApsideDetector(orekit_ss._state)
        at_apoapsis = ApsideDetector(orekit_ss._state).withHandler(StopOnDecreasing())

        # Build the maneuvers
        imp_A = ImpulseManeuver(at_periapsis, attitude_provider, dVa_vec, ISP)
        imp_B = ImpulseManeuver(at_apoapsis, attitude_provider, dVb_vec, ISP)
        burns = (imp_A, imp_B)

        return cls(burns, attitude_provider)
