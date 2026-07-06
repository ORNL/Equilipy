"""Optional MPI helpers for running equilipy batches on HPC clusters.

This module is an add-on: mpi4py is not a base dependency of equilipy.
Install it with the ``hpc`` extra:

    pip install equilipy[hpc]

Design (see _ongoing_development/elegant_minimizer.md, Section 10): batch
equilibrium is embarrassingly parallel over conditions and the Fortran solver
state is per-process, so cross-node scaling is a thin driver layer on top of
the single-node ``equilib_batch`` loky pool.  The recommended hybrid topology
is one MPI rank per node with loky processes inside each node:

    srun --ntasks-per-node=1 --cpus-per-task=<cores> python driver.py

Never call MPI from inside loky workers; MPI stays in the rank-level driver.
"""

from __future__ import annotations

import os

try:
    from mpi4py import MPI
except ImportError as exc:  # pragma: no cover - exercised only without mpi4py
    raise ImportError(
        "equilipy.hpc requires mpi4py, which is not installed. "
        "Install the HPC add-on with: pip install equilipy[hpc]"
    ) from exc

COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
SIZE = COMM.Get_size()


def worker_count() -> int:
    """Return the CPU count available to THIS rank, trusting the scheduler.

    loky's own autodetection sees every core on the node, not the SLURM
    allocation, so relying on it oversubscribes shared nodes.  Precedence:
    SLURM_CPUS_PER_TASK, then the process CPU affinity mask, then
    ``os.cpu_count()``.
    """
    env = os.environ.get("SLURM_CPUS_PER_TASK")
    if env:
        return max(1, int(env))
    if hasattr(os, "sched_getaffinity"):
        try:
            return max(1, len(os.sched_getaffinity(0)))
        except OSError:
            pass
    return os.cpu_count() or 1


def split_condition(condition: dict, rank: int = None, size: int = None) -> dict:
    """Return this rank's contiguous chunk of a condition dictionary.

    Conditions are split into ``size`` near-equal contiguous chunks so that
    gathering per-rank results in rank order preserves the input row order.
    A rank beyond the number of conditions receives an empty chunk.
    """
    if rank is None:
        rank = RANK
    if size is None:
        size = SIZE

    n_total = len(next(iter(condition.values())))
    base, extra = divmod(n_total, size)
    counts = [base + (1 if r < extra else 0) for r in range(size)]
    start = sum(counts[:rank])
    stop = start + counts[rank]
    return {key: values[start:stop] for key, values in condition.items()}
