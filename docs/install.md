# Installation

`equilipy` supports Python 3.10-3.14.

## 1. Install with pip

```bash
pip install equilipy
```

Binary wheels are published for Linux (x86_64, aarch64, and i686; glibc
and musl), Windows (AMD64), and macOS (Apple Silicon) with the Fortran
runtime bundled, so no compiler is required. On other platforms,
including Intel Macs, [install from source](#2-install-from-source).

Optional add-ons are installed as pip *extras*. For the PySide6 desktop GUI
(adds the `equilipy.gui` command):

```bash
pip install 'equilipy[gui]'
```

:::{note}
On minimal Linux systems (servers, cloud images) the GUI can abort with
*"xcb-cursor0 is needed to load the Qt xcb platform plugin"* — Qt needs a
few X11 client libraries that PySide6 does not bundle. Fix it once per
environment: obtain the libraries (e.g.
`conda install -c conda-forge xcb-util-cursor`, or the distro packages),
then run

```bash
equilipy.gui --setup-linux-libs
```

which copies them into PySide6 so a plain `equilipy.gui` works from then
on — no `LD_LIBRARY_PATH` needed.
:::

For the mpi4py helpers on multi-node clusters:

```bash
pip install 'equilipy[hpc]'
```

:::{note}
For the `hpc` extra on a cluster, load the site's MPI module **before**
installing so `mpi4py` compiles against the same MPI that `srun`/`mpirun`
launches with. See [HPC](scripting/hpc) for details.
:::

## 2. Install from source

Building from source — for development, or on platforms without a prebuilt
wheel — requires a Fortran compiler.

With `conda` — Linux:

```bash
conda install -c conda-forge gfortran_linux-64
```

macOS (Intel):

```bash
conda install -c conda-forge gfortran_osx-64
```

macOS (Apple Silicon):

```bash
conda install -c conda-forge gfortran_osx-arm64
```

Windows:

```bash
conda install -c conda-forge fortran-compiler
```

Or with a system package manager — Ubuntu / Debian:

```bash
sudo apt-get install gfortran
```

macOS (Homebrew):

```bash
brew install gcc
```

On Windows without conda, install
[MinGW-w64](https://github.com/niXman/mingw-builds-binaries/releases), copy it
to `C:\mingw\`, and add `C:\mingw\bin` to the `Path` environment variable.

Then build and install the wheel:

```bash
git clone https://github.com/ORNL/Equilipy.git
cd Equilipy
pip install build meson ninja
python -m build --wheel --outdir ./wheelhouse
pip install wheelhouse/equilipy-*.whl
```

On macOS with Homebrew, point the build at GNU `gcc` instead of the Apple
clang `gcc` before building:

```bash
export CC=$(brew --prefix)/bin/gcc-15   # adjust to your installed version
```

:::{note}
When building from source, run the test suite to verify the build:
`pytest tests -q`.
:::
