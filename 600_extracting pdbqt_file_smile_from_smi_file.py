import pandas as pd

# =====================================
# PATHS
# =====================================

SMI_FILE = r"D:\TNBC\28k_smi.csv"
VINA_FILE = r"D:\TNBC\bharati_dock\Result_600_extrcated\unique_best_VINA.csv"
OUTPUT_FILE = r"D:\TNBC\bharati_dock\Result_600_extrcated\600_extracted_smi.csv"

# =====================================
# STEP 1: LOAD SMILES FILE
# =====================================

df_smi = pd.read_csv(SMI_FILE)

# Ensure correct column names
df_smi.columns = ["INDEX", "SMILE"]
df_smi["INDEX"] = pd.to_numeric(df_smi["INDEX"], errors='coerce')

# =====================================
# STEP 2: LOAD VINA FILE
# =====================================

df_vina = pd.read_csv(VINA_FILE)

# Extract numeric ID from Filename column
df_vina["INDEX"] = df_vina["Filename"].str.extract(r'(\d+)').astype(int)

# =====================================
# STEP 3: MERGE DATA
# =====================================

merged_df = pd.merge(df_vina, df_smi, on="INDEX", how="left")

# =====================================
# STEP 4: SAVE RESULT
# =====================================

merged_df.to_csv(OUTPUT_FILE, index=False)

print("\n✅ DONE!")
print(f"Saved at: {OUTPUT_FILE}")
print(f"Total matched: {len(merged_df)}")