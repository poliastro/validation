import matlab
import numpy as np
import pytest
from pytest import approx

from poliastro.core.elements import rv2coe


@pytest.fixture
def earth_elliptic():
    # Taken from poliastro.examples.iss
    r = [8.59072560e2, -4.13720368e3, 5.29556871e3]  # km
    v = [7.37289205, 2.08223573, 4.39999794e-1]  # km / s

    return r, v


@pytest.mark.xfail
def test_rv2coe_iss_exact(engine, earth_elliptic):
    r, v = earth_elliptic

    expected_result = engine.rv2coe(matlab.double(r), matlab.double(v), nargout=11)
    p, _, ecc, inc, omega, argp, nu, *_ = expected_result

    # Load constants to make all computations equivalent
    engine.constastro(nargout=0)
    k = engine.workspace["mu"]

    result = rv2coe(k, np.array(r), np.array(v))

    assert result[0] == p
    assert result[1] == ecc
    assert result[2] == inc
    assert result[3] == omega
    assert result[4] == argp
    assert result[5] == nu


def test_rv2coe_iss_approximate(engine, earth_elliptic):
    r, v = earth_elliptic
    rtol = atol = 1e-13

    expected_result = engine.rv2coe(matlab.double(r), matlab.double(v), nargout=11)
    p, _, ecc, inc, omega, argp, nu, *_ = expected_result

    # Load constants to make all computations equivalent
    engine.constastro(nargout=0)
    k = engine.workspace["mu"]

    result = rv2coe(k, np.array(r), np.array(v))

    assert result[0] == approx(p, rel=rtol, abs=atol)
    assert result[1] == approx(ecc, rel=rtol, abs=atol)
    assert result[2] == approx(inc, rel=rtol, abs=atol)
    assert result[3] == approx(omega, rel=rtol, abs=atol)
    assert result[4] == approx(argp, rel=rtol, abs=atol)
    assert result[5] == approx(nu, rel=rtol, abs=atol)
