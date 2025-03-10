#!/bin/bash
#SBATCH --job-name=MPI1_Tests
#SBATCH --nodes=8
#SBATCH --cpus-per-task=32
#SBATCH --output=bin/mpi1_%j.log
#SBATCH --error=bin/mpi1_%j.err

# Percorso dei dataset
DATASET_PATH_1="/opt/ohpc/pub/SPMcode/LSH/datasets/lsh1GB.dat"
DATASET_PATH_5="/opt/ohpc/pub/SPMcode/LSH/datasets/lsh5GB.dat"
DATASET_PATH_10="/opt/ohpc/pub/SPMcode/LSH/datasets/lsh10GB.dat"

# Percorso dei file di output
OUTPUT_FILE_1="output/mpi1_1GB"
OUTPUT_FILE_5="output/mpi1_5GB"
OUTPUT_FILE_10="output/mpi1_10GB"

# Rimuovi i vecchi file di output
rm -f "$OUTPUT_FILE_1" "$OUTPUT_FILE_5" "$OUTPUT_FILE_10"

# Array delle configurazioni di nodi
NODES=(4 5 6 7 8)

# Array dei dataset e file di output
DATASETS=($DATASET_PATH_1 $DATASET_PATH_5 $DATASET_PATH_10)
OUTPUT_FILES=($OUTPUT_FILE_1 $OUTPUT_FILE_5 $OUTPUT_FILE_10)

# Array delle dimensioni dei file per output leggibile
FILE_SIZES=("1GB" "5GB" "10GB")

# Loop sui dataset
for i in ${!DATASETS[@]}; do
    echo "Starting MPI1 tests on ${FILE_SIZES[$i]} file.." >> "${OUTPUT_FILES[$i]}"
    echo "" >> "${OUTPUT_FILES[$i]}"
    
    # Loop sulle configurazioni dei nodi
    for N in "${NODES[@]}"; do
        echo "Configuration $N nodes:" >> "${OUTPUT_FILES[$i]}"
        srun --mpi=pmix --nodes=$N ./mpi1 "${DATASETS[$i]}" >> "${OUTPUT_FILES[$i]}" 2>&1
        echo "" >> "${OUTPUT_FILES[$i]}"
    done
    
    echo "END OF MPI TESTS for ${FILE_SIZES[$i]}" >> "${OUTPUT_FILES[$i]}"
    echo "" >> "${OUTPUT_FILES[$i]}"

done

echo "Test completati."
