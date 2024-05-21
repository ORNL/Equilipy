---
layout: default
title: Preinstall
nav_order: 4
description: "Before installation"
---

# Before installation
## Single computing nodes (desktop/laptop)
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

## Multiple computing nodes (HPC) on Linux
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
