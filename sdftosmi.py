import os
from rdkit import Chem
import pandas as pd


def sdf_to_smiles(input_folder, output_file):
    """
    Convert all SDF files in a folder to SMILES and save them to a CSV file.

    Parameters:
        input_folder (str): Path to the folder containing SDF files.
        output_file (str): Path to save the output CSV file.
    """
    data = []

    for filename in os.listdir(input_folder):
        if filename.endswith(".sdf"):
            sdf_path = os.path.join(input_folder, filename)

            # Read the SDF file
            suppl = Chem.SDMolSupplier(sdf_path)

            for mol in suppl:
                if mol is not None:
                    smiles = Chem.MolToSmiles(mol)
                    data.append([filename, smiles])

    # Create a DataFrame and save to CSV
    df = pd.DataFrame(data, columns=["Filename", "SMILES"])
    df.to_csv(output_file, index=False)

    print(f"Conversion completed. SMILES saved to {output_file}")


# Example usage
input_folder = 'D:/PROJECT MSC BIOINFO/TNBC/natural drugs/New drugs from databases/impactdrugs'  # Replace with your folder path
output_file = "D:/PROJECT MSC BIOINFO/TNBC/natural drugs/New drugs from databases/imtactsmi.csv"  # Output CSV file path
sdf_to_smiles(input_folder, output_file)
