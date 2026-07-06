"""Input condition transfer helpers for the Fortran backend."""

from __future__ import annotations

import numpy as np

import equilipy.equilifort as fort
import equilipy.variables as var

HEAT_CAPACITY_TEMPERATURE_DIFF = 1e-2


def input_condition(unit, composition, include_heat_capacity=True):
    """Assign input variables in Fortran.

    Description
    ===========
    This function assgins values for input variables in Fortran


    Revisions
    =========

     Date            Programmer      Description of change
     ----            ----------      ---------------------
     11/19/2020      S.Y. KWON       Original code


    Variables
    =========
    units   : A list of units for temperature, pressure, mass e.g.['K','atm','moles']
              Temperature units, 'K'/'C'/'F'/'R'
              Pressure units, 'atm'/'psi'/'bar'/'Pa'/'kPa'
              Mass units, 'grams'/'kilograms'/'pounds'/'g'/'kg'/'lbs'/
                    'mass fraction'/'weight fraction'/'wt%'/'wt.%/
                        'moles'/'mol'/'atoms'/'gram-atoms'/'gram-moles'/
                        'mole fraction'/'atom fraction'/'mol%'/'at%'
    el      : an integer array of atomic number
    comp    : a composition matrix shape of (L,n+2), where L and n are the
              number of compositions and elements each.
    """
    if include_heat_capacity:
        fort.modulethermoio.dtemperaturediff = HEAT_CAPACITY_TEMPERATURE_DIFF
    else:
        fort.modulethermoio.dtemperaturediff = 0.00

    # Assign units
    fort.modulethermoio.cinputunittemperature = "{:<15}".format(unit[0])  # Temperature
    fort.modulethermoio.cinputunitpressure = "{:<15}".format(unit[1])  # Pressure
    fort.modulethermoio.cinputunitmass = "{:<15}".format(unit[2])  # Mass

    # Assign values for temperature and pressure
    fort.modulethermoio.dtemperature = composition[0]
    fort.modulethermoio.dpressure = composition[1]
    composition = composition[2:]

    # Rearanging the composition based on the order given by the database
    if len(composition) > 1:
        composition = composition[var.iElementSysIndex]

    # Assign values for mass
    fort.modulethermoio.delementmass = np.zeros(119)
    fort.modulethermoio.idependentelementinput = np.zeros(119, dtype=np.int32)
    for i in range(len(var.iElementSys)):
        fort.modulethermoio.delementmass[int(var.iElementSys[i])] = composition[i]
    for i in getattr(var, "iPseudoComponentDependentElementSys", []):
        fort.modulethermoio.idependentelementinput[int(i)] = 1

    # Assign value for print option in Fortran
    fort.modulethermoio.iprintresultsmode = 2
    return None
