#!/usr/bin/env python3
"""Run an Equilipy GUI-saved calculation module."""

import os
import equilipy as eq
import polars as pl


fpath = os.path.dirname(os.path.abspath(__file__))

# Module: Equilibrium#1

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
PhaseSelection = ['LIQUID',
 'FCC_A1',
 'FCC_A1#2',
 'HCP_A3',
 'HCP_A3#2',
 'BCC_A2',
 'BCC_A2#2',
 'BCT_A5',
 'BCT_A5#2',
 'DIAMOND_A4',
 'DIAMOND_A4#2',
 'CBCC_A12',
 'CBCC_A12#2',
 'CUB_A13',
 'CUB_A13#2',
 'EPS',
 'EPS#2',
 'SIGMA',
 'SIGMA#2',
 'T_PHASE',
 'T_PHASE#2',
 'GAMMA',
 'GAMMA#2',
 'GAMMA#3',
 'ALCU_THETA',
 'ALCU_THETA#2',
 'ALMG_GAMMA',
 'ALMG_GAMMA#2',
 'EPSILON_CU15SI4',
 'EPSILON_CU15SI4#2',
 'ETA_CU19SI6',
 'ETA_CU19SI6#2',
 'GAMMA_CU56SI11',
 'GAMMA_CU56SI11#2',
 "'AlMg'_ALMG_ZETA(s)",
 "'Al5Mg4'_ALMG_EPSILON(s)",
 "'Al8Mg5'_ALMG_BETA(s)",
 "'Mg2Al11'_MG2ZN11(s)",
 'Mg2Si_MG2SI(s)',
 "'Mg3Si'_MG19SI6(s)",
 "'Mg4Si'_MG15SI4(s)",
 "'Mg5Si'_MG56SI11(s)"]
NPT = {'T': 1600.0, 'P': 1.0, 'Al': 93.6, 'Mg': 0.4000000000000001, 'Si': 6.0}

# Step 3: Calculate equilibrium
res = eq.equilib_single(
    DB,
    NPT,
    unit=UnitIn,
    phases=PhaseSelection,
)

# Step 4: Save and print result
result_filename = 'Equilibrium_1.eqres'
result_path = os.path.join(fpath, result_filename)
if hasattr(res, 'to_bundle'):
    eq.save_result(res, result_path)
res_table = res.to_dict() if hasattr(res, 'to_dict') else res
print(pl.DataFrame(res_table) if isinstance(res_table, dict) else res_table)
