# IsosurfaceCpp

A C++ library consisting of different isosurface extraction algorithms.
Its only dependency is the mathematics library [glm](https://glm.g-truc.net/0.9.9/).

An example of how to use IsosurfaceCpp as a submodule can be found on https://github.com/chrismile/LineVis.
The library can also be compiled and installed as a standalone package using CMake.


### Python Module

IsosurfaceCpp also supports building a Python module. The dependencies for the Python module can be installed via `pip`:

```sh
pip install pybind11 numpy setuptools wheel
```

After all dependencies have been installed, the Python module can be built with `pip` using the command specified below.
It seems like `--no-build-isolation` is only necessary when installing in a Python venv, not a conda environment.

```sh
pip install --no-build-isolation .
```

An example for the usage of the Python module can be found in `pymodule/example.py`.
