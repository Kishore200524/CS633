#!/bin/bash
#SBATCH --job-name=a2_P32_n240
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:10:00
#SBATCH --output=slurm_outputs/P32_n240_akki_%j.out
#SBATCH --error=slurm_outputs/P32_n240_akki_%j.err

# Configuration: P=32, px=4, py=4, pz=2, ppn=32, nx=ny=nz=240

module load compiler/oneapi-2024/mpi

EXE=./akki_assg2
D=7; PPN=32; PX=4; PY=4; PZ=2; NX=240; NY=240; NZ=240; T=5; SEED=1000; F=2; ISO=500

for run in 1 2 3 4 5; do
    echo "=== Run $run ==="
    mpirun -np 32 $EXE $D $PPN $PX $PY $PZ $NX $NY $NZ $T $SEED $F $ISO
done
