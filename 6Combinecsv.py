import pandas as pd
import os

# === Step 1: Combine all CSV files ===
input_directory = r"E:\TNBC\admet\Rushi\Download"
combined_output_file = r"E:\TNBC\admet\Rushi\Download\combined_.csv"
filtered_output_file = r"E:\TNBC\admet\Rushi\Download\filtered_molecules.csv"

# List all CSV files in the directory
csv_files = [f for f in os.listdir(input_directory) if f.endswith('.csv')]

# Initialize an empty DataFrame
combined_df = pd.DataFrame()

# Loop through all CSV files and concatenate them
for file in csv_files:
    file_path = os.path.join(input_directory, file)
    df = pd.read_csv(file_path)
    combined_df = pd.concat([combined_df, df], ignore_index=True)

# Save the combined DataFrame to a single CSV file
combined_df.to_csv(combined_output_file, index=False)
print(f"Combined {len(csv_files)} CSV files into '{combined_output_file}'")

# === Step 2: Apply filtering criteria ===
# Load the combined CSV file
df = pd.read_csv(combined_output_file)

# Apply the filtering conditions
filtered_df = df[
    (df['QED'] > 0.61) &
    (df['PAINS'] == "['-']") &
    (df['Lipinski'] == 0) &
    (df['GoldenTriangle'] == 0) &
    (df['caco2'] > -5.15) &
    (df['BBB'] < 0.31) &
    (df['hERG'] < 0.31) &
    (df['DILI'] < 0.31) &
    (df['Ames'] < 0.31) &
    (df['H-HT'] < 0.6) &
    (df['FAF-Drugs4 Rule'] == "['-']")
]

# Save the filtered results
filtered_df.to_csv(filtered_output_file, index=False)
print(f"Filtering complete. Results saved to '{filtered_output_file}'")
