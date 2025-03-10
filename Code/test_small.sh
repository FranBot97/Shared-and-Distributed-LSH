#!/bin/bash
#SBATCH --job-name=Small_test
#SBATCH --nodes=8
#SBATCH --cpus-per-task=16
#SBATCH --output=bin/mpi_%j.log
#SBATCH --error=bin/mpi_%j.err

# Percorso dei dataset
DATASET_PATH="/opt/ohpc/pub/SPMcode/LSH/datasets/lsh1GB.dat"
OUTPUT_FILE="output/small.out"

#Command example
#srun --nodes=1 --ntasks=1 ./fastflow /opt/ohpc/pub/SPMcode/LSH/datasets/lsh1GB.dat 7 6
#srun --nodes=1 --ntasks=1 ./sequential "/opt/ohpc/pub/SPMcode/LSH/datasets/lsh1GB.dat"


rm -f "$OUTPUT_FILE"

echo "Starting all the small tests on 1GB file" >> "$OUTPUT_FILE" 2>&1
echo "" >> "$OUTPUT_FILE"

echo "Sequential 1GB" >> "$OUTPUT_FILE"
srun --nodes=1 --ntasks=1 ./sequential "${DATASET_PATH}" >> "$OUTPUT_FILE" 2>&1
echo "" >> "$OUTPUT_FILE"

echo "Fastflow (7,6) 1GB" >> "$OUTPUT_FILE"
srun --nodes=1 --ntasks=1 ./fastflow "${DATASET_PATH}" 7 6 >> "$OUTPUT_FILE" 2>&1
echo "" >> "$OUTPUT_FILE"	

echo "MPI1 8 nodes" >> "$OUTPUT_FILE"
srun --mpi=pmix --nodes=8 ./mpi1 "${DATASET_PATH}" >> "$OUTPUT_FILE" 2>&1
echo "" >> "$OUTPUT_FILE"

echo "MPI2 8 nodes" >> "$OUTPUT_FILE"
srun --mpi=pmix --nodes=8 ./mpi2 "${DATASET_PATH}" >> "$OUTPUT_FILE" 2>&1
echo "" >> "$OUTPUT_FILE"


echo "END OF SMALL TESTS" >> "$OUTPUT_FILE"
echo "Test completati."
