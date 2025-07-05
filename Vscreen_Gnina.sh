#!/bin/bash

# Paths
LIGAND_DIR="./pdb_files"
RECEPTOR=$(ls ./*.pdbqt | head -n 1)  # Get the first .pdbqt file in the current directory
CONFIG="./config.txt"
OUTPUT_DIR="./Output"

# Check if a receptor file is found
if [ -z "$RECEPTOR" ]; then
    echo "Error: No .pdbqt receptor file found in the current directory."
    exit 1
fi

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Iterate through all ligand files in the ligand directory
for LIGAND in "$LIGAND_DIR"/*.pdbqt; do
    # Get base name of ligand file (without extension)
    BASE_NAME=$(basename "$LIGAND" .pdbqt)
    
    echo -e "\n\nDocking ligand ---> $BASE_NAME\n"
    
    # Create output directory for each ligand
    LIGAND_OUTPUT_DIR="$OUTPUT_DIR/$BASE_NAME"
    mkdir -p "$LIGAND_OUTPUT_DIR"
    
    # Run Gnina using Docker
    sudo docker run --rm \
        -v "$(pwd)":/data \
        -v "$LIGAND_DIR":/data/Inhibitors \
        -v "$OUTPUT_DIR":/data/Output \
        gnina/gnina gnina \
        -r "/data/$(basename "$RECEPTOR")" \
        -l "/data/Inhibitors/$(basename "$LIGAND")" \
        -o "/data/Output/$BASE_NAME/${BASE_NAME}_out.pdbqt" \
        --config "/data/$(basename "$CONFIG")"
done
