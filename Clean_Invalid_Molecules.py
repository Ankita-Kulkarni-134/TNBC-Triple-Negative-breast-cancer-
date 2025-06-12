from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

# Function to validate SMILES
def validate_smiles(input_file):
    valid_smiles = []
    invalid_smiles = []

    with open(input_file, "r") as infile:
        for line in infile:
            smiles = line.strip()
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                try:
                    # Using MorganGenerator instead of deprecated method
                    morgan_gen = AllChem.GetMorganGenerator(2)  # 2 is the radius
                    morgan_fingerprint = morgan_gen.GetFingerprint(mol)
                    valid_smiles.append(smiles)
                except Exception as e:
                    print(f"Error processing SMILES '{smiles}': {e}")
                    invalid_smiles.append(smiles)
            else:
                invalid_smiles.append(smiles)

    return valid_smiles, invalid_smiles

# Path to the original SMILES file
original_file_path = r"C:\Users\admin\Desktop\TNBC\admet\Rushi\Combined_unique.smi"

# Validate SMILES in the original file
valid_smiles, invalid_smiles = validate_smiles(original_file_path)

# Write valid SMILES to a new file in SMILES format
with open(r"C:\Users\admin\Desktop\TNBC\admet\Rushi\Combined_unique_valid.smi", "w") as outfile:
    outfile.write("\n".join(valid_smiles))

# Write invalid SMILES to a new file in SMILES format
with open(r"C:\Users\admin\Desktop\TNBC\admet\Rushi\Combined_unique_invalid.smi", "w") as outfile:
    outfile.write("\n".join(invalid_smiles))

# Now, filter to keep only valid SMILES
filtered_valid_smiles = valid_smiles  # No need to read the file again

# Write the filtered valid SMILES to a new file in SMILES format
with open(r"C:\Users\admin\Desktop\TNBC\admet\Rushi\filtered_valid_smiles.smi", "w") as outfile:
    outfile.write("\n".join(filtered_valid_smiles))
