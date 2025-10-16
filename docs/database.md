---
layout: default
title: Database
nav_enabled: true
nav_order: 5
---

# Thermochemical Database

## Database format
`equilipy` uses the data format (`.dat`) provided by ChemSage (also known as ChemApp). Detailed information about the database format can be found in  [ChemApp_online-manual](https://gtt-technologies.de/software/chemapp/documentation/online-manual/). This data format is compatible with commercial programs like [FactSage](https://www.factsage.com/) and non-commercial programs such as [PyCalphad](https://pycalphad.org/docs/latest/) and [Thermochimica](https://github.com/ORNL-CEES/thermochimica).

The current version of `equilipy` does not support `.tdb` database format. However, the conversion from `.tdb` to `.dat` is available in [FactSage7.3](https://www.factsage.com/).

{: .warning }
The ChemSage `.dat` format was recently changed in [FactSage8.0+](https://www.factsage.com/).  The new `.dat` format is not yet compatible with `equilipy`.

## Example database
`equilipy` provides an example thermochemical database for the Al-Cu-Mg-Si quaternary system. Thermodynamic assessments of partinent unary, binary, and ternary systems were taken from [COST507](https://materialsdata.nist.gov/handle/11256/618). 


