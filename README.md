# Equilipy
Equilipy is an open-source python package that offers multicomponent multiphase equilibrium calculations based on the CALPHAD (CALculation of PHAse Diagram) approach. With a set of Gibbs energy description (Thermochemical database) and input conditions (Composition, temperature, pressure), equilibrium phase configureation, amount, composition, and thermochemical properties can be obtained. Equilipy uses the Gibbs energy descriptions furnished by THERMOCHIMICA with the modified Gibbs energy minimization algorithm initially proposed by de Capitani, C. and Brown, T.H. (1987).

## Dependencies
|Dependency | Version  | Required | Libraries |
|---------- | -------  |--------  |-------    |
|Fortran    | -        | Yes      | -
|Python     | 3.9+     | Yes      | numpy, wheel, meson, ninja

## Before installation
### Single computing nodes (desktop/laptop)
Equilipy requires a Fortran compiler in the local environment.
To install gfortran using `conda`,
for **Linux**:
```
conda install -c conda-forge gfortran_linux-64
```
for **MacOS**:
```
conda install -c conda-forge gfortran_osx-64
```
for **Windows**:
```
conda install -c conda-forge fortran-compiler
```

Alternatively, gfortran can be install for **Ubuntu** and **Debian**,
```
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install gfortran 
sudo apt-get install libopenmpi-dev
```
for **MacOS**,
```
brew install gcc open-mpi
```
To install gfortran on **Windows**,
1. Download the latest [MinGW-w64](https://github.com/niXman/mingw-builds-binaries/releases) and unzip.
2. Copy the unzipped folder to C-drive and rename the folder/directory as **mingw** in C-drive `C:\mingw\`
3. Click the **Windows** button and type "environment variables" to access **Edit the system environment variables**.
4. Click the **Environment Variables** at the bottom right corner
5. Click **Path** in System variables dialog to display **Edit environment variable** window
6. Click **New** and add `C:\mingw\bin` to the path

Equilipy also requires Python version 3.9 and above. The Fortran backend needs to be compiled through the f2py module in `numpy` which requires `meson` and `ninja`. The `wheel` library is used for packaging. These can both be installed through pip.

### Multiple computing nodes (HPC) on Linux
Equilipy uses `mpi4py` to interface with MPI tools. 
To install OpenMPI, mpi4py, and gfortran without using **sudo** privilage, we recommand install gfortran, OpenMPI, and mpi4py using `conda`:
```
conda install -c conda-forge gfortran_linux-64 openmpi mpi4py
```

Alternatively, users with **sudo** privilage may install without `conda`:
for **Debian-based** (Debian, Ubuntu, Mint, etc..)
```
sudo apt-get install gfortran 
sudo apt-get install libopenmpi-dev
```
for **CentOS** (Red Hat Enterprise Linux, CentOS,Fedora, openSUSE),
```
sudo yum install gcc-gfortran 
sudo yum install openmpi openmpi-devel
```
for **Fedora** and **Red Hat Enterprise Linux**
```
sudo dnf install gcc-fortran
sudo dnf install openmpi openmpi-devel
```

## Installation

### Pypi
Installation using `pip` is available for Equilipy.
```
pip install equilipy
```

For HPC environment, use `equilipy-hpc`
```
pip install equilipy-hpc
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