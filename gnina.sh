#!/bin/bash

# Navigate to the docking directory
cd /mnt/d/PROJECTS/TNBC/multi_dock_458/multi_dock_458

# Make sure output folder exists
mkdir -p results

# Loop through each ligand (.pdbqt) file
for ligand in Ligands/*.pdbqt; do
    # Extract the base file name (without extension)
    base=$(basename "$ligand" .pdbqt)

    # Run GNINA docking with your config file
    sudo docker run -v $(pwd):/data gnina/gnina gnina \
    -r /data/6q4g.pdbqt \
    -l "/data/Ligands/${base}.pdbqt" \
    -o "/data/results/${base}_out.sdf" \
    --config /data/config.txt

    echo "Docking done for: $base"
done

echo " All 458 ligands have been docked successfully!"
