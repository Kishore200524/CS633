#!/bin/bash
#SBATCH --job-name=assign1_p32
#SBATCH --nodes=2              
#SBATCH --ntasks=32            
#SBATCH --ntasks-per-node=16   
#SBATCH --time=00:10:00
#SBATCH --output=output1_p32.txt

module load compiler/oneapi-2024/mpi

echo "========================================"
echo "Starting runs for P=32"
echo "========================================"

# Run for Small M
echo "----------------------------------------"
echo "Running for M = 262144 (Small)"
echo "----------------------------------------"
for i in {1..5}
do
    mpirun -np 32 ./src 262144 2 4 10 1000
done

# Run for Large M
echo ""
echo "----------------------------------------"
echo "Running for M = 1048576 (Large)"
echo "----------------------------------------"
for i in {1..5}
do
    mpirun -np 32 ./src 1048576 2 4 10 1000
done