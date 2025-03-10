#!/bin/bash
#SBATCH --job-name=Sequential_Test
#SBATCH --nodes=3  # Richiedi 3 nodi per l'esecuzione parallela
#SBATCH --cpus-per-task=32
#SBATCH --output=bin/sequential_%j.log
#SBATCH --error=bin/sequential_%j.err

# Percorso dei dataset
DATASET_PATH_1="/opt/ohpc/pub/SPMcode/LSH/datasets/lsh1GB.dat"
DATASET_PATH_5="/opt/ohpc/pub/SPMcode/LSH/datasets/lsh5GB.dat"
DATASET_PATH_10="/opt/ohpc/pub/SPMcode/LSH/datasets/lsh10GB.dat"

OUTPUT_FILE_1="output/sequential_1GB.out"
OUTPUT_FILE_5="output/sequential_5GB.out"
OUTPUT_FILE_10="output/sequential_10GB.out"

rm -f "$OUTPUT_FILE_1"
rm -f "$OUTPUT_FILE_5"
rm -f "$OUTPUT_FILE_10"

# Esegui i test in parallelo su nodi diversi
echo "1GB:" >> "$OUTPUT_FILE_1"
srun --nodes=1 --exclusive ./sequential "${DATASET_PATH_1}" >> "$OUTPUT_FILE_1" 2>&1 &

echo "5GB:" >> "$OUTPUT_FILE_5"
srun --nodes=1 --exclusive ./sequential  "${DATASET_PATH_5}" >> "$OUTPUT_FILE_5" 2>&1 &

echo "10GB:" >> "$OUTPUT_FILE_10"
srun --nodes=1 --exclusive ./sequential  "${DATASET_PATH_10}" >> "$OUTPUT_FILE_10" 2>&1 &

# Attendi la fine di tutti i processi in background
wait