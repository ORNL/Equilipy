# Equilipy
Equilipy is an open-source python package that offers multicomponent multiphase equilibrium calculations based on the CALPHAD (CALculation of PHAse Diagram) approach. With a set of Gibbs energy description (Thermochemical database) and input conditions (Composition, temperature, pressure), equilibrium phase configureation, amount, composition, and thermochemical properties can be obtained. Equilipy uses the Gibbs energy descriptions furnished by THERMOCHIMICA with the modified Gibbs energy minimization algorithm initially proposed by de Capitani, C. and Brown, T.H. (1987).

## Dependencies
|Dependency | Version  | Required | Libraries |
|---------- | -------  |--------  |-------    |
|Fortran    | -        | Yes      | -
|Python     | 3.9+     | Yes      | numpy, wheel, meson, ninja

## Before installation
Because Equilipy uses Fortran subroutines a Fortran compiler needs to be available in the local envirionment. To install GNU Fortran, following command may be used for Ubuntu and Debian,
```
sudo apt-get install gfortran
```
for MacOS,
```
brew install gcc
```
for Windows, check out [MinGW-w64](https://www.mingw-w64.org/downloads/#mingw-builds).

Equilipy also requires Python version 3.9 and above. The Fortran backend needs to be compiled through the f2py module in `numpy` which requires `meson` and `ninja`. The `wheel` library is used for packaging. These can both be installed through pip:
```
pip install numpy wheel meson ninja
```

## Installation
Here is the procedure to install Equilipy:
1. Clone the repository
```
git clone https://github.com/ORNL/Equilipy.git
```
2. Create wheels
```
python setup.py bdist_wheel --dist-dir=./wheelhouse
```
Note that on macOS it may be necessary to explicitly use GNU `gcc` instead of the Apple clang `gcc`. For example, with a brew installation, export `CC` before calling `setup.py`:
```
export CC=/usr/local/Cellar/gcc/*/bin/gcc-*
```
3. Install from .whl file
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
