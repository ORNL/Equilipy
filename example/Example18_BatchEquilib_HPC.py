#!/usr/bin/env python3
"""HPC hybrid (MPI x loky) variant of Example03_BatchEquilib.

Requires the HPC add-on: pip install equilipy[hpc]

Topology: one MPI rank per node; each rank solves a contiguous chunk of the
input conditions with the ordinary single-node ``eq.equilib_batch`` (loky
process pool), then rank 0 gathers the per-rank tables and writes one CSV.

Run on a SLURM cluster (see submit_Example18_HPC.sh):

    srun --ntasks-per-node=1 --cpus-per-task=<cores> \
        python Example18_BatchEquilib_HPC.py
"""
import os
import time
from datetime import timedelta

import polars as pl

import equilipy as eq
from equilipy.hpc import COMM, RANK, SIZE, split_condition, worker_count

system = "AlCuMgSi"

# Step 1: Parse database. Each rank parses locally; solver state is per-process.
fpath = os.path.dirname(os.path.abspath(__file__))
path = os.path.dirname(fpath)
DB = eq.read_dat(os.path.join(path, "database", "AlCuMgSi_ORNL_FS83.dat"))

# Step 2: Read input conditions and take this rank's contiguous chunk.
NPT = pl.read_excel(os.path.join(fpath, "Input_ACMS.xlsx"), sheet_name=system).to_dict()
chunk = split_condition(NPT)
n_local = len(chunk["T"])

# Step 3: Solve this rank's chunk with the ordinary single-node batch.
starttime = time.time()
res = eq.equilib_batch(
    DB,
    chunk,
    unit=["K", "atm", "moles"],
    n_cpu=worker_count(),
    progress=(RANK == 0),
)
duration = time.time() - starttime
print(
    f"[rank {RANK}/{SIZE}] {n_local} conditions, "
    f"n_cpu={worker_count()}, time={timedelta(seconds=duration)}"
)

# Step 4: Gather per-rank tables on rank 0 (rank order preserves row order)
# and write one CSV. Chunks can expose different stable-phase columns, so the
# concat must be diagonal.
tables = COMM.gather(res.to_dict() if n_local else None, root=0)
if RANK == 0:
    df = pl.concat([pl.DataFrame(t) for t in tables if t], how="diagonal")
    print(df)
    df.write_csv(f"{fpath}/Result_Ex03_{system}_HPC.csv")
