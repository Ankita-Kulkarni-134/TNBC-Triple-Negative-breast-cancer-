import pandas as pd
from rdkit import Chem
from tqdm import tqdm
import sascorer

# =====================================
# PATHS
# =====================================

INPUT_FILE = r"D:\TNBC\bharati_dock\Result_600_extrcated\600_extracted_smi.csv"
OUTPUT_FILE = r"D:\TNBC\bharati_dock\Result_600_extrcated\600_with_SA.csv"

# =====================================
# STEP 1: LOAD DATA
# =====================================

df = pd.read_csv(INPUT_FILE)

# Check column name
# print(df.columns)

SMI_COL = "SMILE"   # based on your file

# =====================================
# STEP 2: CALCULATE SA SCORE
# =====================================

sa_scores = []

print("\n🔹 Calculating Synthetic Accessibility...\n")

for smi in tqdm(df[SMI_COL]):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            score = sascorer.calculateScore(mol)
        else:
            score = None
    except:
        score = None

    sa_scores.append(score)

df["SA_Score"] = sa_scores

# =====================================
# STEP 3: SAVE OUTPUT
# =====================================

df.to_csv(OUTPUT_FILE, index=False)

print("\n✅ DONE!")
print(f"Saved at: {OUTPUT_FILE}")