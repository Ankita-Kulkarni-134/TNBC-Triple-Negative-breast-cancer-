import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# Disable RDKit warnings
RDLogger.DisableLog('rdApp.*')

# ================================
# INPUT FILE
# ================================
INPUT_SMI = r"D:\TNBC\28k_smi.smi"

# ================================
# OUTPUT DIRECTORIES
# ================================
PDB_DIR = r"D:\TNBC\PDB_files"
PDBQT_DIR = r"D:\TNBC\PDBQT_files"

# ================================
# ✅ OPENBABEL PATH (CONFIRMED WORKING)
# ================================
OBABEL_PATH = r"C:\Program Files\OpenBabel-2.4.1\obabel.exe"

# CPU cores
NUM_CORES = 10

# Create folders
os.makedirs(PDB_DIR, exist_ok=True)
os.makedirs(PDBQT_DIR, exist_ok=True)


def process_molecule(data):
    idx, smi = data

    try:
        pdb_path = os.path.join(PDB_DIR, f"mol_{idx}.pdb")
        pdbqt_path = os.path.join(PDBQT_DIR, f"mol_{idx}.pdbqt")

        # ✅ Skip already processed molecules (resume support)
        if os.path.exists(pdbqt_path):
            return None

        # Convert SMILES → Mol
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return f"Invalid SMILES at {idx}"

        mol = Chem.AddHs(mol)

        # 3D embedding
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True

        result = AllChem.EmbedMolecule(mol, params)
        if result != 0:
            return f"Embedding failed at {idx}"

        # Save PDB
        Chem.MolToPDBFile(mol, pdb_path)

        # 🔥 Convert PDB → PDBQT (FIXED for Open Babel 2.4.1)
        conversion = subprocess.run(
            [
                OBABEL_PATH,
                "-ipdb", pdb_path,
                "-opdbqt",
                "-O", pdbqt_path
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        # Check for conversion errors
        if conversion.returncode != 0:
            return f"OpenBabel failed at {idx}: {conversion.stderr.decode(errors='ignore')}"

        # Ensure file was created
        if not os.path.exists(pdbqt_path):
            return f"PDBQT not created at {idx}"

        return None

    except Exception as e:
        return f"Error at {idx}: {str(e)}"


def main():
    # Check OpenBabel path
    if not os.path.exists(OBABEL_PATH):
        print("❌ OpenBabel not found at given path!")
        print("👉 Fix OBABEL_PATH before running.")
        return

    # Read SMILES file
    with open(INPUT_SMI, "r") as f:
        smiles_list = [line.strip().split()[0] for line in f if line.strip()]

    print(f"Total molecules found: {len(smiles_list)}")
    print(f"Using {NUM_CORES} CPU cores...\n")

    failed = []

    # Parallel processing
    with ProcessPoolExecutor(max_workers=NUM_CORES) as executor:
        results = list(tqdm(
            executor.map(process_molecule, enumerate(smiles_list, 1)),
            total=len(smiles_list)
        ))

    # Collect failures
    for r in results:
        if r is not None:
            failed.append(r)

    print("\n✅ Conversion completed.")
    print(f"❌ Failed molecules: {len(failed)}")

    # Save failure log
    if failed:
        log_path = r"D:\TNBC\failed_smiles.txt"
        with open(log_path, "w") as f:
            for item in failed:
                f.write(item + "\n")

        print(f"📝 Failure log saved at: {log_path}")


if __name__ == "__main__":
    main()