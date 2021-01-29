import os

from orekit.pyhelpers import download_orekit_data_curdir

if not os.path.exists("orekit-data.zip"):
    try:
        download_orekit_data_curdir()
    except Exception as e:
        raise RuntimeError("Orekit data was not successfully downloaded") from e
