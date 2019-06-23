import matlab
import numpy as np


def test_dummy(engine):
    assert engine.isreal(matlab.double([1]))
    assert engine.isinteger(matlab.int8([1]))


def test_path_is_correct(engine):
    # Keep nargout=0, otherwise we get
    # matlab.engine.MatlabExecutionError: Attempt to execute SCRIPT constastro as a function
    engine.constastro(nargout=0)

    assert engine.workspace["twopi"] == 2 * np.pi
