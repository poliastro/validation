# Orekit validation framework

The following directory holds a validation framework against the [Orekit Python
wrapper](https://gitlab.orekit.org/orekit-labs/python-wrapper).


## Users installation

Preferred way of installation, due to the nature of Orekit, is through `conda`.
Thus, an `environment.yml` file is provided, listing all required dependencies
for this project. Build a new environment and install required packages by
running:


```bash
conda env create -f environment.yml
``` 

Finally, activate the environment by running:

```bash
conda activate poliastro-validation
```

where `poliastro-validation` is the name of the previously created conda
environment.


## Developers installation

If you require to use the validation utilities for testing some poliastro source
code modification, then you should install poliastro in development in the conda
environment. Therefore, start by creating the conda environment by running:

```bash
conda create -n poliastro-validation
``` 

and activate it using:

```bash
conda activate poliastro-validation
```

Install the required dependencies via:

```
conda install -c conda-forge flit 
conda install -c conda-forge orekit
conda install -c anaconda pytest
```

then clone the poliastro repository and install the library in development mode
by following [the official
guide](https://docs.poliastro.space/en/stable/contributing.html#development-environment)

## How to run the tests

The `orekit/tests/` directory holds all necessary configuration files and unit
tests to be read during the validation process. Therefore, simply run:

```bash
pytest tests/
```

Previous command will start a `pytest` session who's configuration is hosted in
`orekit/tests/pytest.ini`.


## A detailed installation of Orekit wrapper

Although the [Orekit](https://gitlab.orekit.org/orekit/orekit) original library
is developed under Java programming language, a Python wrapper is provided to
easy interface with all the different routines and classes of the package.

There exist two different ways to [install the
wrapper](https://gitlab.orekit.org/orekit-labs/python-wrapper/-/wikis/installation):
automatic or manual. 

1. The automatic installation requires from a Anaconda installation, since
   Orekit packages are built in the conda-forge system. Then, simply run:

   ```bash
   conda install -c conda-forge orekit
   ```
   
   to install the wrapper and all additional requirements such us openJDK
   together with necessary path variables.

2. Regarding the manual installation process, a detailed explanation is
   available within the
   [INSTALL.txt](https://gitlab.orekit.org/orekit-labs/python-wrapper/blob/master/INSTALL.txt)
   file in the official Orekit wrapper repository.


## About orekit-data.zip

In the official installation guide of the Orekit wrapper, it is stated that
physical data (such us epehemeris) needs to be [downlaoded from official
sources](https://www.orekit.org/download.html).

This file is loaded by the `setup_orekit_curdir()`, which by default will look
for a ZIP file under the name `orekit-data.zip`.

