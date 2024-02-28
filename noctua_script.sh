#!/bin/bash

#SBATCH -t 30:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -J test
#SBATCH -A hpc-prf-flucurhi
#SBATCH --array=0-2

module load numlib/GSL/2.7-GCC-11.3.0
module load toolchain/intel/2022.00
module load compiler/GCC/11.2.0
module load devel/CMake/3.18.4-GCCcore-10.2.0
module load mpi/OpenMPI/4.0.5-GCC-10.2.0
module load lang
module load Anaconda3

conda activate myenv

INPUT_DIR="input_energy_momentum_tensors"
FILES=("$INPUT_DIR"/*.dat)

CURRENT_FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"

if [ ! -f "$CURRENT_FILE" ]; then
    echo "File not found: $CURRENT_FILE"
    exit 1
fi

./ExecuteEBE.sh "$CURRENT_FILE"
