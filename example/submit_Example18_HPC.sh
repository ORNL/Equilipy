#!/bin/bash
#SBATCH --job-name=eq_ex03
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1        # one MPI rank per node (hybrid: loky inside the node)
#SBATCH --cpus-per-task=32         # cores handed to each rank's loky pool
#SBATCH --time=01:00:00
#SBATCH --output=eq_ex03_%j.out

# --- One-time environment setup (login node, not in this script) -----------
# mpi4py on PyPI is a source build: it compiles against whatever MPI `mpicc`
# points to.  Load the site's MPI module BEFORE installing, so mpi4py links
# against the same MPI that srun launches with:
#
#   module load openmpi            # or the site's MPI module
#   source activate equilipy
#   pip install equilipy[hpc]
#
# Alternative: if the site provides a prebuilt mpi4py module (built against
# the interconnect-tuned MPI), prefer it and install only the base package:
#
#   module load python mpi4py
#   pip install equilipy
# ----------------------------------------------------------------------------

# Runtime: load the SAME modules mpi4py was built against.
# module load openmpi
# source activate equilipy

# The Fortran extension is single-threaded per process; keep BLAS/OpenMP from
# oversubscribing under the loky workers.
export OMP_NUM_THREADS=1

srun python Example18_BatchEquilib_HPC.py
