# Equilipy
Equilipy is an open-source python package that offers multicomponent multiphase equilibrium calculations based on the CALPHAD (CALculation of PHAse Diagram) approach. With a set of Gibbs energy description (Thermochemical database) and input conditions (Composition, temperature, pressure), equilibrium phase configureation, amount, composition, and thermochemical properties can be obtained. Equilipy uses the Gibbs energy descriptions furnished by THERMOCHIMICA with the modified Gibbs energy minimization algorithm initially proposed by de Capitani, C. and Brown, T.H. (1987).

## Dependencies
|Dependency | Version  | Required | Libraries |
|---------- | -------  |--------  |-------    |
|Fortran    | -        | Yes      | -
|Python     | 3.9+     | Yes      | numpy, wheel, meson, ninja

## Before installation
Equilipy requires a Fortran compiler and MPI in the local envirionment. To install GNU Fortran and OpenMPI in Linux using `conda`

```
conda install -c conda-forge gfortran_linux-64 openmpi mpi4py
```
or in MacOS
```
conda install -c conda-forge gfortran_osx-64 openmpi mpi4py
```

Alternatively, following command may be used for Ubuntu and Debian,
```
sudo apt-get install gfortran libopenmpi-dev
```
for MacOS,
```
brew install gcc open-mpi
```
for Windows, MS-MPI is pre-installed. For GNU Fortran, check out [MinGW-w64](https://www.mingw-w64.org/downloads/#mingw-builds).

Equilipy also requires Python version 3.9 and above. The Fortran backend needs to be compiled through the f2py module in `numpy` which requires `meson` and `ninja`. The `wheel` library is used for packaging. These can both be installed through pip.

## Installation

### Pypi
Installation using `pip` is available for Equilipy. Please ensure OpenMPI (or MPICH) is preinstalled:
```
pip install equilipy
```

### Install from the source
Here is the procedure to install Equilipy from the source:
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

## Features and example
The following features are currently available.
- Single condition equilibrium calculations
- Batch equilibrium calculations
- Scheil-Gulliver solidification
- Phase selection

For details, check out the example directory.
Running the first example is as simple as:
```
cd example/
python Example01_SingleEquilib.py
```

## Additional note
Examples in Equilipy uses `polars` dataframe for fast data processing. If you are using old CPUs, install
```
pip install polars-lts-cpu
```
If you are using large dataset (> 4billion), install 
```
pip install polars-u64-idx
```
For details, check out [polars](https://docs.pola.rs/).