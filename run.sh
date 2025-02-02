#!/bin/bash

# Request resources:
#SBATCH -n 16       # 1 CPU core
#SBATCH --mem=8G      # 1 GB RAM
#SBATCH --time=0-48:0:0  # 1 hour (days-hours:minutes:seconds)
#SBATCH -p shared
#SBATCH --array=1-1

# Make python available:
module load python
module load gcc
module load openmpi
module load netcdf

pip install numpy

# Commands to be run:

echo "Task number $SLURM_ARRAY_TASK_ID"

python run_smartpressure.py ${SLURM_ARRAY_TASK_ID}

mpirun -np 16 ./bin/mf3d ${SLURM_ARRAY_TASK_ID}

