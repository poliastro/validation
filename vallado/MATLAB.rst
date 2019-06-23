Using the MATLAB Python bridge
==============================

Installation
------------

To install the MATLAB Python engine using pipenv:

1. Give write permissions in
   :code:`/usr/local/MATLAB/R2019a/extern/engines/python/dist`.
2. Run :code:`pipenv install -e /usr/local/MATLAB/R2019a/extern/engines/python/`.
3. Create a
   :code:`/usr/local/MATLAB/R2019a/extern/engines/python/dist/matlab/engine/_arch.txt`
   with the contents below.
4. Profit!

Contents of :code:`_arch.txt`::

  glnxa64
  /usr/local/MATLAB/R2019a/extern/engines/python/../../../bin/glnxa64
  /usr/local/MATLAB/R2019a/extern/engines/python/dist/matlab/engine/glnxa64

Caveats:

* Alternatively you can mangle :code:`PYTHONPATH` to find the shared library directly
* Notice that the official method uses distutils
  and therefore makes it impossible to uninstall
* Did not test our method with pip but it should work
* Do not try to use non editable mode, it fails (perhaps because it tries to build a wheel?)
* The method described in
  https://www.scivision.dev/matlab-engine-callable-from-python-how-to-install-and-setup/
  to avoid giving write permissions in the MATLAB directory
  does not create the required `_arch.txt` file
* You might have to fix libssl as described in
  https://github.com/imatlab/imatlab/issues/3#issuecomment-289574029
  when trying to run pytest
  (perhaps because https://github.com/pyenv/pyenv/issues/1112)

Usage
-----

Basic usage example::

  import matlab.engine
  eng = matlab.engine.start_matlab("-nodesktop -nojvm")

  x = matlab.double([1.0, 0.5])
  y = eng.asin(x)

  print(y)

  eng.quit()

Returning several output arguments::

  import matlab.engine
  eng = matlab.engine.start_matlab("-nodesktop -nojvm")

  t = eng.gcd(100, 80, nargout=3)
  print(t)

  eng.quit()

Catching errors and redirecting output::

  import io
  import matlab.engine
  eng = matlab.engine.start_matlab("-nodesktop -nojvm")

  out = io.StringIO()
  err = io.StringIO()

  try:
      eng.dec2base(2 ** 60, 16, stdout=out, stderr=err)
  except matlab.engine.MatlabExecutionError:
      pass

  print(out.getvalue())
  print(err.getvalue())

  eng.quit()

Other
-----

* To put Python variables into MATLAB workspace see
  https://www.mathworks.com/help/matlab/apiref/matlab.engine.matlabengine-class.html#buki0wz
