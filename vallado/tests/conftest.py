import os

import matlab.engine
import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
MATLAB_SOURCES_DIR = os.path.join(TESTS_DIR, os.pardir, "matlab")


# TODO: Use wider scope? Tests would be more fragile
@pytest.fixture()
def engine():
    eng = matlab.engine.start_matlab("-nodesktop -nojvm")
    eng.addpath(MATLAB_SOURCES_DIR)
    yield eng
    eng.quit()
