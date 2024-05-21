---
layout: default
title: Home
nav_order: 1
description: "Quick description"
---

# Computing multicomponent multiphase equilibria

Equilipy is an open-source python package that offers multicomponent multiphase equilibrium calculations based on the [CALPHAD][CALPHAD method] (CALculation of PHAse Diagram) approach. With a set of Gibbs energy description (Thermochemical database) and input conditions (Composition, temperature, pressure), equilibrium phase configureation, amount, composition, and thermochemical properties can be obtained. Equilipy uses the Gibbs energy descriptions furnished by [THERMOCHIMICA][Thermochimica] with the modified Gibbs energy minimization algorithm initially proposed by de Capitani, C. and Brown, T.H.[^1]

[View it on GitHub][equilipy]{: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }

---

## Getting started

For details, check out [Features and example][features].
Running the first example is as simple as:

## Quick Installation

Installation using `pip` is available for Equilipy.
```
pip install equilipy
```

{: .note }
Examples in Equilipy uses `polars` dataframe for fast data processing. If you are using old CPUs, use ```pip install polars-lts-cpu```. If you are using large dataset (> 4billion), use 
```pip install polars-u64-idx```. For details, check out [Polars][polars].

For detailed information, please see [Install][install].

{: .warning }
> Equilipy requires a Fortran compiler in the local environment. To install a Fortran compiler, please follow [Preinstall][preinstall].

## About the project

Euqilipy is &copy; 2024-{{ "now" | date: "%Y" }} by [U.S. Department of Energy](https://doi.org/10.11578/dc.20240312.4).

### License

Euqilipy is distributed by an [BSD 3-Clause License](https://github.com/ORNL/Equilipy/blob/main/LICENSE).

### Code of Conduct

We as contributors and maintainers pledge to making participation in our project and
our community a poisitive experiences for everyone.

[View our Code of Conduct](https://github.com/ORNL/Equilipy/blob/main/CODE_OF_CONDUCT.md) on our GitHub repository.

----
[^1]: de Capitani, C. and Brown, T.H., *Geochimica et Cosmochimica Acta* 51, (1987): 2639-2652.

[equilipy]: https://github.com/ORNL/Equilipy
[CALPHAD method]: https://www-sciencedirect-com.ornl.idm.oclc.org/journal/calphad/about/aims-and-scope
[Thermochimica]: https://github.com/ORNL-CEES/thermochimica
[preinstall]: https://ornl.github.com/Equilipy/blob/main/docs/preinstal.mdl
[install]: https://ornl.github.com/Equilipy/blob/main/docs/install.md
[features]: https://ornl.github.com/Equilipy/blob/main/docs/features.md
[polars]: https://docs.pola.rs/