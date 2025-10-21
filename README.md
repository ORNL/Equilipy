[![PyPI - Version](https://img.shields.io/pypi/v/equilipy)](https://pypi.org/project/equilipy/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/equilipy)](https://pypistats.org/packages/equilipy)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.06875/status.svg)](https://doi.org/10.21105/joss.06875)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17408401.svg)](https://doi.org/10.5281/zenodo.17408401)

# Equilipy
Equilipy is an open-source python package that offers multicomponent multiphase equilibrium calculations based on the CALPHAD (CALculation of PHAse Diagram) approach. With a set of Gibbs energy description (Thermochemical database) and input conditions (Composition, temperature, pressure), equilibrium phase configureation, amount, composition, and thermochemical properties can be obtained. Equilipy uses the Gibbs energy descriptions furnished by THERMOCHIMICA with the modified Gibbs energy minimization algorithm initially proposed by de Capitani, C. and Brown, T.H. (1987).

Check out [documentation](https://ornl.github.io/Equilipy/) for further description.

## Dependencies
|Dependency | Version  | Required | Libraries |
|---------- | -------  |--------  |-------    |
|Fortran    | -        | Yes      | -
|Python     | 3.9+     | Yes      | numpy, wheel, meson, ninja


## Installation

Installation using `pip` is available for Equilipy.
```
pip install equilipy
```

## Features and example
The following features are currently available.
- Single condition equilibrium calculations
- Batch equilibrium calculations
- Scheil-Gulliver solidification
- Phase selection
Hehehe

For details, check out the example directory and [Features and Examples](https://ornl.github.io/Equilipy/features.html)

## Citing Equilipy

If you use Equilipy in your work, please cite the following [paper](CITATION.bib).
In addition, cite the current release or version used from
[Zenodo](https://doi.org/10.5281/zenodo.17408401).

## Contributing
We encourage you to contribute to Equilipy. Please see [contributing guidelines](CONTRIBUTING.md).

## Additional note
Examples in Equilipy uses `polars` dataframe for fast data processing. In particular, example 3 requires `fastexcel` as the optional dependancy in `polars`. 
Install `fastexcel` via 
```
pip install fastexcel
```
Additionally, if you are using large dataset (> 4billion), install 
```
pip install polars-u64-idx
```
If you are using old CPUs, install
```
pip install polars-lts-cpu
```
For details, check out [polars dependencies](https://docs.pola.rs/api/python/stable/reference/api/polars.show_versions.html).