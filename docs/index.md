# Equilipy for computing multicomponent-multiphase equilibria

`equilipy` is an open-source Python package for multicomponent-multiphase
equilibrium calculations based on the [CALPHAD][CALPHAD method]
(CALculation of PHAse Diagram) approach. Given a thermochemical database and
input conditions (composition, pressure, and temperature), it computes the
equilibrium phase configuration, amounts, compositions, and thermochemical
properties. `equilipy` uses the Gibbs energy descriptions furnished by
[THERMOCHIMICA][Thermochimica] with the modified Gibbs energy minimization
algorithm inspired by Eriksson[^1] and de Capitani and Brown.[^2]

```{toctree}
:maxdepth: 2
:caption: User guide

install
gui/index
scripting/index
```

% Developer documentation (architecture.md, developer/) is kept in the repo
% but excluded from the published build via exclude_patterns in conf.py.
% To publish it again, remove the exclusions and restore the toctree:
%
% ```{toctree}
% :maxdepth: 2
% :hidden:
%
% architecture
% developer/index
% ```

## Quick start

### Quick installation

Core package:

```bash
pip install equilipy
```

Desktop GUI (adds the `equilipy.gui` command):

```bash
pip install 'equilipy[gui]'
```

MPI helpers for clusters:

```bash
pip install 'equilipy[hpc]'
```

% Verified against the published 0.2.1 wheels (all bundle libgfortran via
% delocate/auditwheel/delvewheel), but keep this hidden until the no-compiler
% install is re-verified with the v0.3.0 release wheels:
%
% :::{note}
% Wheels on PyPI ship precompiled binaries with the Fortran runtime bundled,
% so `pip install` needs no compiler. Only building from source requires a
% Fortran compiler; see [Installation](install).
% :::

### Database availability

Equilipy reads two thermochemical database formats:

| Format | Reader | Notes |
|---|---|---|
| ChemSage/FactSage `.dat` | `eq.read_dat(path)` | Legacy (7.3) and modern (8.3+) dialects are auto-detected, including validated MQM/SUBG/SUBQ solution-model canaries. |
| Thermo-Calc style `.tdb` | `eq.read_tdb(path)` | CEF solution models, magnetic contributions, and order/disorder (`DIS_PART`) are supported. |

An example Al-Cu-Mg-Si database (assessments from
[COST507](https://www.metallurgy.nist.gov/reports/nistir6927.pdf)) ships in the
repository's `database/` folder. For larger databases such as [Al–Co–Cr–Fe–Mn–Ni–C](https://doi.org/10.1016/j.calphad.2023.102644), please check out publications in [CALPHAD journal](https://www.sciencedirect.com/journal/calphad).

### Example

```python
import equilipy as eq

DB = eq.read_dat("database/AlCuMgSi_ORNL_FS83.dat")

NPT = {"T": 700, "P": 1, "Al": 0.06, "Cu": 0.42, "Si": 0.52}
res = eq.equilib_single(DB, NPT)

print(res.stable_phases.names)   # stable phase names
print(res.G, res.H, res.S)       # system properties, J and J/K
```

More runnable scripts are in the repository's [example folder][examples]:

```bash
python example/Example01_SingleEquilib.py
```

Continue with [Python scripting](scripting/index) for batch calculations,
solidification, and results handling, or [GUI](gui/index) for the desktop
application.

:::{note}
Examples use the [Polars][polars] dataframe library, installed by default.
On old CPUs use `pip install polars-lts-cpu`.
:::

## About the project

Equilipy is copyright 2024-present by [U.S. Department of Energy](https://doi.org/10.11578/dc.20240312.4).

### License

Equilipy is distributed under a [BSD 3-Clause License](https://github.com/ORNL/Equilipy/blob/main/LICENSE).

### Code of Conduct

We as contributors and maintainers pledge to make participation in our project
and community a positive experience for everyone.

[View our Code of Conduct](https://github.com/ORNL/Equilipy/blob/main/CODE_OF_CONDUCT.md) on our GitHub repository.

[^1]: Eriksson, G., *Acta Chemica Scandinavica* 25, (1971): 2651-2658.
[^2]: de Capitani, C. and Brown, T.H., *Geochimica et Cosmochimica Acta* 51, (1987): 2639-2652.

[equilipy]: https://github.com/ORNL/Equilipy
[CALPHAD method]: https://calphad.org
[Thermochimica]: https://github.com/ORNL-CEES/thermochimica
[polars]: https://docs.pola.rs/
[examples]: https://github.com/ORNL/Equilipy/blob/main/example
