---
layout: default
title: Install
nav_enabled: true
nav_order: 3
---

# Installation

## PIP
Installation using `pip` is available for `equilipy`. 
To install `equilipy` in your **desktop/laptop**:
```
pip install equilipy
```

{: .warning }
`equilipy` requires a Fortran compiler in the local environment. To install a Fortran compiler, please follow [Preinstall][preinstall].

{: .note }
The test suite in the `./test/*` directory is integrated into the build process for .whl files using GitHub Actions. No further testing is necessary when installing via pip.

For **HPC environment**, use `equilipy-hpc` instead of `equilipy`
```
pip install equilipy-hpc
```

{: .warning }
`equilipy-hpc` uses `mpi4py` to interface with MPI tools. MPI tools such as OpenMPI or MPICH must be preinstalled. To install To install OpenMPI and mpi4py, check out [Preinstall][preinstall].

## Install from the source
[GitHub repository][equilipy] provides source codes to build `equilipy`.
To install `equilipy` from the source:
1. Clone the repository
```
git clone https://github.com/ORNL/Equilipy.git
```
2. Install dependant packages using pip.
```
pip install numpy wheel meson ninja
```
3. Create wheels
```
python setup.py bdist_wheel --dist-dir=./wheelhouse
```
Note that on macOS it may be necessary to explicitly use GNU `gcc` instead of the Apple clang `gcc`. For example, with a brew installation, export `CC` before calling `setup.py`:
```
export CC=/usr/local/Cellar/gcc/*/bin/gcc-*
```
4. Install from .whl file
```
pip install wheelhouse/equilipy-*.whl
```

{: .note }
When building equilipy from the source, consider running the test suite located in the `./test/*` directory. This test suite ensures the reliability of the core functions in equilipy.

[preinstall]: https://github.com/ORNL/Equilipy/blob/main/docs/preinstall.md
[equilipy]: https://github.com/ORNL/Equilipy
