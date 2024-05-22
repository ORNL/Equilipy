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

For details, check out the example directory and [Features and Examples](https://ornl.github.io/Equilipy/features.html)


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