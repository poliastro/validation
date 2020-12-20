import os
import sys

import pytest

# Get orekit validation tests/ and find the src/ one
TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
OREKIT_SRC_DIR = os.path.join(TESTS_DIR, os.pardir, "src")

# Append src/ to path to access custom orekit modules
sys.path.append(OREKIT_SRC_DIR)
