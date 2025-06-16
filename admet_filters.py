import pandas as pd
import os
from difflib import get_close_matches
from tqdm import tqdm

# === Input and Output folder paths ===
input_folder = r"D:\TNBC\ADMET_CSV"  # Folder with 2000 CSV files
output_folder = r"D:\TNBC\Admet_filtered_outputs"  # Output folder for filtered CSVs

# === Ensure output folder exists ===
os.makedirs(output_folder, exist_ok=True)

# === Define numeric columns for filtering ===
numeric_columns = ['QED', 'Lipinski', 'GoldenTriangle', 'BBB',
                   'Ames', 'hERG', 'H-HT', 'DILI']

# === Summary list ===
summary_list = []

# === Process each CSV file with progress bar ===
for filename in tqdm(os.listdir(input_folder), desc="Processing CSV files"):
    if filename.endswith(".csv"):
        input_path = os.path.join(input_folder, filename)
        output_path = os.path.join(output_folder, f"filtered_{filename}")

        try:
            # Load the CSV
            df = pd.read_csv(input_path)

            # Normalize column names
            df.columns = df.columns.str.strip().str.replace('\u200b', '').str.replace('\n', ' ').str.replace('\r', '').str.strip()

            # Fuzzy match the 'FAF-Drugs4 rule' column
            target_col = 'FAF-Drugs4 rule'
            matches = get_close_matches(target_col, df.columns, n=3, cutoff=0.5)
            if matches:
                df.rename(columns={matches[0]: target_col}, inplace=True)
            else:
                summary_list.append((filename, "SKIPPED - FAF-Drugs4 rule not found", 0))
                continue

            # Ensure numeric columns are valid
            df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors='coerce')
            df_cleaned = df.dropna(subset=numeric_columns)

            # Normalize string columns
            df_cleaned['PAINS'] = df_cleaned['PAINS'].astype(str).str.strip().str.lower()
            df_cleaned['FAF-Drugs4 rule'] = df_cleaned['FAF-Drugs4 rule'].astype(str).str.strip().str.lower()

            # DEBUG: Print unique values of PAINS and FAF-Drugs4 rule
            print(f"\n{filename} - Unique PAINS values: {df_cleaned['PAINS'].unique()}")
            print(f"{filename} - Unique FAF values: {df_cleaned['FAF-Drugs4 rule'].unique()}")

            # Define relaxed numeric filter conditions
            conditions = {
                'QED': df_cleaned['QED'] > 0.5,
                'Lipinski': df_cleaned['Lipinski'] <= 1,
                'GoldenTriangle': df_cleaned['GoldenTriangle'] <= 1,
                'BBB': df_cleaned['BBB'] < 0.5,
                'Ames': df_cleaned['Ames'] < 0.5,
                'hERG': df_cleaned['hERG'] < 0.5,
                'H-HT': df_cleaned['H-HT'] < 0.7,
                'DILI': df_cleaned['DILI'] < 0.5,
            }

            # Score each molecule
            df_cleaned['FilterScore'] = sum(cond.astype(int) for cond in conditions.values())

            # DEBUG: Print score distribution
            print(f"{filename} - Filter score distribution:\n", df_cleaned['FilterScore'].value_counts())

            # Filter based on relaxed score only (no PAINS/FAF rule for now)
            filtered_df = df_cleaned[df_cleaned['FilterScore'] >= 4].drop(columns=['FilterScore'])

            # DEBUG: Print number of rows passing
            print(f"{filename} - Rows passing relaxed filter: {len(filtered_df)}")

            # Save the filtered results
            filtered_df.to_csv(output_path, index=False)

            # Log result
            summary_list.append((filename, "PROCESSED", len(filtered_df)))

        except Exception as e:
            summary_list.append((filename, f"ERROR - {e}", 0))

# === Save summary report ===
summary_df = pd.DataFrame(summary_list, columns=["Filename", "Status", "Filtered_Row_Count"])
summary_df.to_csv(os.path.join(output_folder, "filtering_summary.csv"), index=False)

print("\n✅ Filtering complete.")
print(f"🔍 Summary saved to: {os.path.join(output_folder, 'filtering_summary.csv')}")
