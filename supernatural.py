import os
import requests
import pandas as pd

# Input CSV file with compound IDs
input_file = r"/mnt/data/indication_Breast_cancer.csv"
output_file = r"/mnt/data/output_smiles.smi"

# Load compound IDs from the input CSV
data = pd.read_csv(input_file)
compound_ids = data.iloc[:, 0]  # Assuming the first column contains IDs

# Base URL for the SuperNatural3.0 database
base_url = "https://bioinf-applied.charite.de/supernatural_3/search/"

# Open the output file for writing
with open(output_file, "w") as outfile:
    for compound_id in compound_ids:
        try:
            # Query the database for SMILES
            response = requests.get(f"{base_url}?id={compound_id}")
            response.raise_for_status()

            # Parse the response for SMILES (assuming JSON format with 'smiles' key)
            smiles = response.json().get('smiles', None)

            if smiles:
                outfile.write(f"{compound_id}\t{smiles}\n")
                print(f"Extracted SMILES for: {compound_id}")
            else:
                print(f"No SMILES found for: {compound_id}")
        except requests.exceptions.RequestException as e:
            print(f"Failed to fetch data for {compound_id}: {e}")
