#!/bin/bash
#SBATCH --job-name=a2_P96_n120
#SBATCH --nodes=2
#SBATCH --ntasks=96
#SBATCH --ntasks-per-node=48
#SBATCH --time=00:10:00
#SBATCH --output=slurm_outputs/P96_n120_akki_%j.out
#SBATCH --error=slurm_outputs/P96_n120_akki_%j.err

# Configuration: P=96, px=6, py=4, pz=4, ppn=48, nx=ny=nz=120

module load compiler/oneapi-2024/mpi

EXE=./akki_assg2
D=7; PPN=48; PX=6; PY=4; PZ=4; NX=120; NY=120; NZ=120; T=5; SEED=1000; F=2; ISO=500

for run in 1 2 3 4 5; do
    echo "=== Run $run ==="
    mpirun -np 96 $EXE $D $PPN $PX $PY $PZ $NX $NY $NZ $T $SEED $F $ISO
done
