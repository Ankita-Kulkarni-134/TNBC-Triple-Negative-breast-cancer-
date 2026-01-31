import pandas as pd
import os
from difflib import get_close_matches
from tqdm import tqdm

# ================================
# INPUT & OUTPUT PATHS
# ================================
input_folder = r"D:\PROJECTS\TNBC\sample_admet_out"
output_folder = r"D:\TNBC\sample_admet_out2"

os.makedirs(output_folder, exist_ok=True)

print("Saving filtered files to:", output_folder)

# ================================
# NUMERIC COLUMNS
# ================================
numeric_columns = [
    'QED',
    'Lipinski',
    'GoldenTriangle',
    'caco2',
    'BBB',
    'Ames',
    'DILI'
]

# ================================
# STORAGE FOR ALL FILTERED DATA
# ================================
all_filtered_dfs = []

# ================================
# SUMMARY LOG
# ================================
summary_list = []

# ================================
# PROCESS FILES
# ================================
for filename in tqdm(os.listdir(input_folder), desc="Processing CSV files"):
    if not filename.lower().endswith(".csv"):
        continue

    input_path = os.path.join(input_folder, filename)
    output_path = os.path.join(output_folder, f"filtered_{filename}")

    try:
        # ----------------------------
        # LOAD CSV
        # ----------------------------
        df = pd.read_csv(input_path)

        # ----------------------------
        # CLEAN COLUMN NAMES
        # ----------------------------
        df.columns = (
            df.columns.astype(str)
            .str.replace('\u200b', '')
            .str.replace('\n', ' ')
            .str.replace('\r', '')
            .str.strip()
        )

        # ----------------------------
        # FUZZY MATCH FAF-DRUGS4 RULE
        # ----------------------------
        target_col = 'FAF-Drugs4 rule'
        matches = get_close_matches(target_col, df.columns, n=3, cutoff=0.5)

        if matches:
            df.rename(columns={matches[0]: target_col}, inplace=True)
        else:
            summary_list.append((filename, "SKIPPED - FAF-Drugs4 rule not found", 0))
            continue

        # ----------------------------
        # REQUIRED COLUMNS CHECK
        # ----------------------------
        required_cols = numeric_columns + ['PAINS', target_col]
        missing_cols = [c for c in required_cols if c not in df.columns]

        if missing_cols:
            summary_list.append((filename, f"SKIPPED - Missing columns: {missing_cols}", 0))
            continue

        # ----------------------------
        # NUMERIC CONVERSION
        # ----------------------------
        df[numeric_columns] = df[numeric_columns].apply(
            pd.to_numeric, errors='coerce'
        )

        df_cleaned = df.dropna(subset=numeric_columns)

        if df_cleaned.empty:
            summary_list.append((filename, "SKIPPED - No valid rows after numeric cleaning", 0))
            continue

        # ----------------------------
        # NORMALIZE STRING COLUMNS
        # ----------------------------
        df_cleaned = df_cleaned.copy()

        df_cleaned['PAINS'] = (
            df_cleaned['PAINS']
            .astype(str)
            .str.strip()
            .str.lower()
        )

        df_cleaned[target_col] = (
            df_cleaned[target_col]
            .astype(str)
            .str.strip()
            .str.lower()
        )

        # ----------------------------
        # FILTER CONDITIONS
        # ----------------------------
        conditions = {
            'QED': df_cleaned['QED'] > 0.5,
            'PAINS': df_cleaned['PAINS'] == '[-]',
            'Lipinski': df_cleaned['Lipinski'] == 1,
            'GoldenTriangle': df_cleaned['GoldenTriangle'] == 0,
            'caco2': df_cleaned['caco2'] > -5.15,
            'BBB': df_cleaned['BBB'] < 0.31,
            'Ames': df_cleaned['Ames'] < 0.31,
            'DILI': df_cleaned['DILI'] < 0.31
        }

        # ----------------------------
        # SCORE
        # ----------------------------
        df_cleaned['FilterScore'] = sum(
            cond.astype(int) for cond in conditions.values()
        )

        # ----------------------------
        # RELAXED FILTER
        # ----------------------------
        filtered_df = (
            df_cleaned[df_cleaned['FilterScore'] >= 4]
            .drop(columns=['FilterScore'])
        )

        if filtered_df.empty:
            summary_list.append((filename, "PROCESSED - No rows passed filter", 0))
            continue

        # ----------------------------
        # SAVE INDIVIDUAL FILE
        # ----------------------------
        filtered_df.to_csv(output_path, index=False)

        # ----------------------------
        # STORE FOR GLOBAL CONCAT
        # ----------------------------
        filtered_df['Source_File'] = filename
        all_filtered_dfs.append(filtered_df)

        summary_list.append((filename, "PROCESSED", len(filtered_df)))

    except Exception as e:
        summary_list.append((filename, f"ERROR - {str(e)}", 0))

# ================================
# CONCATENATE ALL FILTERED RESULTS
# ================================
if all_filtered_dfs:
    combined_df = pd.concat(all_filtered_dfs, ignore_index=True)

    combined_output_path = os.path.join(
        output_folder, "ALL_FILTERED_COMBINED.csv"
    )

    combined_df.to_csv(combined_output_path, index=False)

    print(f"\n✅ Combined filtered file saved to:\n{combined_output_path}")
    print(f"🧪 Total molecules after filtering: {len(combined_df)}")
else:
    print("\n⚠️ No filtered data available to combine.")

# ================================
# SAVE SUMMARY
# ================================
summary_df = pd.DataFrame(
    summary_list,
    columns=["Filename", "Status", "Filtered_Row_Count"]
)

summary_path = os.path.join(output_folder, "filtering_summary.csv")
summary_df.to_csv(summary_path, index=False)

print("\n✅ Filtering complete.")
print(f"🔍 Summary saved to: {summary_path}")
