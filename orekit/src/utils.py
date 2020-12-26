""" Utilities for working with Orekit wrapper """

import logging
import os

import requests

import orekit
from orekit.pyhelpers import setup_orekit_curdir

logger = logging.getLogger(__name__)

OREKIT_SRC_DIR = os.path.dirname(os.path.abspath(__file__))
OREKIT_DATA_DIR = os.path.join(OREKIT_SRC_DIR, "data")
OREKIT_DATA_PATH = os.path.join(OREKIT_DATA_DIR, "orekit-data.zip")


def build_data_dir(dirname_path=OREKIT_DATA_DIR):
    """ Builds orekit data directory """

    try:
        os.mkdir(dirname_path)
    except FileExistsError:
        pass


def download_orekit_data(save_path=OREKIT_DATA_PATH, chunk_size=128):
    """ Downloads orekit ZIP data """

    # Official data download link
    url = "https://gitlab.orekit.org/orekit/orekit-data/-/archive/master/orekit-data-master.zip"

    # Start downloading process and inform user when done
    r = requests.get(url, stream=True)
    with open(save_path, "wb") as fd:
        for chunk in r.iter_content(chunk_size=chunk_size):
            fd.write(chunk)
    logger.info("Data downloaded correctly.")


def setup_orekit_env(data_path=OREKIT_DATA_PATH):
    """ Sets up the orekit virtual machine and data directory """

    # Create data directory if necessary
    if not os.path.exists(OREKIT_DATA_DIR):
        build_data_dir(OREKIT_DATA_DIR)

    # Download data ZIP file if required
    if not os.path.exists(OREKIT_DATA_PATH):
        download_orekit_data(OREKIT_DATA_PATH)

    # Setup the virtual machine and data location
    vm = orekit.initVM()
    setup_orekit_curdir(OREKIT_DATA_PATH)

    return vm


vm = setup_orekit_env()
""" Sets up Orekit's virtual machine """
