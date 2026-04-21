#!/usr/bin/env python3

import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import RDLogger

# 🔇 Disable RDKit warnings (important)
RDLogger.DisableLog('rdApp.*')

# ── File paths ─────────────────────────────────────────────
INPUT_FILE = r"D:\TNBC\bharati_dock\Result\all_molecules.smi"
OUTPUT_CSV = r"D:\TNBC\bharati_dock\smi_mw_results.csv"


def process_smi():
    seen = set()
    results = []

    with open(INPUT_FILE, 'r') as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            smiles = parts[0]
            name = parts[1] if len(parts) > 1 else f"mol_{i}"

            # 🔁 Remove duplicates
            if smiles in seen:
                continue
            seen.add(smiles)

            # 🧪 Convert SMILES → molecule
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                results.append((name, smiles, "", "Invalid SMILES"))
                continue

            # ⚖️ Calculate MW
            mw = Descriptors.MolWt(mol)

            results.append((name, smiles, round(mw, 3), ""))

    return results


def main():
    print("Processing SMILES file...")

    results = process_smi()

    # ── Save CSV ─────────────────────────────────────────
    with open(OUTPUT_CSV, 'w') as f:
        f.write("Sr_No,Name,SMILES,Molecular_Weight,Error\n")

        for i, (name, smiles, mw, err) in enumerate(results, 1):
            f.write(f"{i},{name},{smiles},{mw},{err}\n")

    print(f"\n✅ Done! CSV saved at:\n{OUTPUT_CSV}")
    print(f"Total molecules processed: {len(results)}")


if __name__ == "__main__":
    main()