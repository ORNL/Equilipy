# HPC

Batch equilibrium is embarrassingly parallel over conditions, so scaling to
a cluster is a thin MPI driver on top of the ordinary `equilib_batch`: one
MPI rank per node, and inside each node the usual process pool.

## Install

Load the site's MPI module first, then install with the `hpc` extra so
`mpi4py` compiles against the same MPI that `srun` launches with:

```bash
module load openmpi          # or the site's MPI module
pip install 'equilipy[hpc]'
```

If the site provides a prebuilt `mpi4py` module, prefer it and install the
base package only (`module load python mpi4py`, then `pip install equilipy`).

## Hybrid driver

`equilipy.hpc` provides the rank-level helpers:

```python
import polars as pl

import equilipy as eq
from equilipy.hpc import COMM, RANK, SIZE, split_condition, worker_count

DB = eq.read_dat("AlCuMgSi_ORNL_FS83.dat")          # each rank parses locally

NPT = pl.read_excel("Input_ACMS.xlsx", sheet_name="AlCuMgSi").to_dict()  # fastexcel engine, installed with equilipy
chunk = split_condition(NPT)                         # this rank's rows

res = eq.equilib_batch(
    DB,
    chunk,
    n_cpu=worker_count(),                            # cores of THIS rank
    progress=(RANK == 0),                            # one progress bar
)

tables = COMM.gather(res.to_dict(), root=0)          # rank order preserved
if RANK == 0:
    df = pl.concat([pl.DataFrame(t) for t in tables if t], how="diagonal")
    df.write_csv("results.csv")
```

- `worker_count()` trusts the scheduler (`SLURM_CPUS_PER_TASK`, then CPU
  affinity) instead of the raw node core count, so shared nodes are not
  oversubscribed.
- `split_condition` slices contiguous near-equal chunks, so gathering in
  rank order preserves the input row order.
- The concat is *diagonal* because different chunks can expose different
  stable-phase columns.

A complete runnable script is
[`example/Example18_BatchEquilib_HPC.py`](https://github.com/ORNL/Equilipy/blob/main/example/Example18_BatchEquilib_HPC.py).

## SLURM submission

```bash
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1        # one MPI rank per node
#SBATCH --cpus-per-task=32         # cores handed to each rank's pool

module load openmpi
export OMP_NUM_THREADS=1           # solver is single-threaded per process

srun python Example18_BatchEquilib_HPC.py
```

A template is provided as
[`example/submit_Example18_HPC.sh`](https://github.com/ORNL/Equilipy/blob/main/example/submit_Example18_HPC.sh).

:::{note}
For plain parameter sweeps, a SLURM job array splitting the condition table
by index gives the same scaling with no MPI dependency. Use MPI when you
want a single output file or in-run coordination.
:::
