#!/usr/bin/env python3
"""Run an Equilipy GUI-saved calculation module."""

import os
import equilipy as eq
import polars as pl


fpath = os.path.dirname(os.path.abspath(__file__))

# Module: Solidification#1

# Step 1: Parse database
database_path = 'databases/AlCuMgSi_ORNL_FS83.dat'
if not database_path:
    raise RuntimeError('No database path was saved with this module.')
if not os.path.isabs(database_path):
    database_path = os.path.join(fpath, database_path)
if database_path.lower().endswith('.tdb'):
    DB = eq.read_tdb(database_path, strict=False, auto_correct=False)
else:
    DB = eq.read_dat(database_path)

# Step 2: Set input data
UnitIn = ['C', 'atm', 'wt%']
PhaseSelection = None
NPT = {'T': 2000.0, 'P': 1.0, 'Al': 92.3, 'Mg': 5.5, 'Si': 2.2}
TargetPhase = 'LIQUID'

# Step 3: Calculate solidification
res = eq.scheil_cooling(
    TargetPhase,
    DB,
    NPT,
    delta_T=5.0,
    unit=UnitIn,
    phases=PhaseSelection,
    start_from_liquidus=True,
)

# Step 4: Save and print result
result_filename = 'Solidification_1.eqres'
result_path = os.path.join(fpath, result_filename)
if hasattr(res, 'to_bundle'):
    eq.save_result(res, result_path)
res_table = res.to_dict() if hasattr(res, 'to_dict') else res
print(pl.DataFrame(res_table) if isinstance(res_table, dict) else res_table)
