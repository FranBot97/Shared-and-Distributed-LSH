#!/bin/bash
#SBATCH --job-name=Fastflow_Test
#SBATCH --nodes=3  # Richiede 3 nodi
#SBATCH --cpus-per-task=32
#SBATCH --output=bin/fastflow_%j.log
#SBATCH --error=bin/fastflow_%j.err

# Percorso dei dataset
declare -A DATASET_PATHS=(
    [1]="/opt/ohpc/pub/SPMcode/LSH/datasets/lsh1GB.dat"
    [5]="/opt/ohpc/pub/SPMcode/LSH/datasets/lsh5GB.dat"
    [10]="/opt/ohpc/pub/SPMcode/LSH/datasets/lsh10GB.dat"
)

CONFIGURATIONS=(1 2 4 5 6 8 10 14)

# Funzione per eseguire i test in sequenza su un nodo specifico
run_tests() {
    local dataset_size=$1
    local dataset_path=${DATASET_PATHS[$dataset_size]}
    local output_file="output/fastflow_${dataset_size}GB.out"

    rm -f "$output_file"
    echo "Starting FastFlow tests on ${dataset_size}GB file.." >> "$output_file"
    echo "" >> "$output_file"

    for config in "${CONFIGURATIONS[@]}"; do
        echo "Configuration (${config},${config}):" >> "$output_file"
        srun --nodes=1 --ntasks=1 --wait=0 ./fastflow "$dataset_path" "$config" "$config" >> "$output_file" 2>&1
        echo "" >> "$output_file"
    done

    echo "END OF TESTS for ${dataset_size}GB" >> "$output_file"
}

# Esegui i test su 3 nodi in parallelo, ma ogni nodo li esegue in sequenza
run_tests 1 &  # Nodo 1 ? 1GB
run_tests 5 &  # Nodo 2 ? 5GB
run_tests 10 & # Nodo 3 ? 10GB

wait  # Aspetta che tutti i nodi abbiano terminato

echo "All tests completed."
