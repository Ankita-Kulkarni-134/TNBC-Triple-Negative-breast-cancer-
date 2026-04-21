import pandas as pd
import shutil
import os

# ── paths ──────────────────────────────────────────────────────────────────
XLSX_PATH   = r"D:\TNBC\bharati_dock\less_than_600.xlsx"
SOURCE_DIRS = [
    r"D:\TNBC\bharati_dock\6q4gbatch_1\src\Output",
    r"D:\TNBC\bharati_dock\6q4gbatch_2\src\Output",
    r"D:\TNBC\bharati_dock\batch3\Output"
]
DEST_DIR    = r"D:\TNBC\bharati_dock\Result_600_extrcated"

# ── load filenames from first column ───────────────────────────────────────
df = pd.read_excel(XLSX_PATH)
first_col = df.columns[0]          # "filename" column
filenames = df[first_col].dropna().astype(str).str.strip().tolist()

# Ensure every name ends with .pdbqt
target_files = set()
for name in filenames:
    if not name.lower().endswith(".pdbqt"):
        name = name + ".pdbqt"
    target_files.add(name)

print(f"Total files to look for: {len(target_files)}")

# ── create destination folder ───────────────────────────────────────────────
os.makedirs(DEST_DIR, exist_ok=True)

# ── search & copy ──────────────────────────────────────────────────────────
found     = []
not_found = list(target_files)   # start assuming none found

for src_dir in SOURCE_DIRS:
    if not os.path.isdir(src_dir):
        print(f"[WARNING] Source folder not found, skipping: {src_dir}")
        continue

    for fname in os.listdir(src_dir):
        if fname in target_files:
            src_path  = os.path.join(src_dir, fname)
            dest_path = os.path.join(DEST_DIR, fname)

            shutil.copy2(src_path, dest_path)
            found.append(fname)
            if fname in not_found:
                not_found.remove(fname)
            print(f"  ✔  Copied: {fname}")

# ── summary ────────────────────────────────────────────────────────────────
print(f"\n{'='*50}")
print(f"✅  Copied   : {len(found)} file(s)")
print(f"❌  Not found: {len(not_found)} file(s)")

if not_found:
    print("\nFiles listed in Excel but NOT found in any source folder:")
    for f in sorted(not_found):
        print(f"   - {f}")