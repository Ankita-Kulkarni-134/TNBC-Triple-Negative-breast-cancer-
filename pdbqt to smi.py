#!/usr/bin/env python3
"""
PDBQT to SMI Converter
========================
Converts all .pdbqt files in the input folder to .smi (SMILES) format.

Input  folder : D:\TNBC\bharati_dock\Result
Output folder : D:\TNBC\bharati_dock\SMI_Output

Requirements:
    pip install openbabel-wheel
    OR install OpenBabel separately and ensure it's in PATH.

Run:
    python pdbqt_to_smi.py
"""

import os
import sys
import subprocess
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
INPUT_FOLDER  = r"D:\TNBC\bharati_dock\Result"
OUTPUT_FOLDER = r"D:\TNBC\bharati_dock\SMI_Output"

# ── Conversion settings ───────────────────────────────────────────────────────
# If True  → one combined .smi file with all molecules
# If False → individual .smi file per .pdbqt file
COMBINED_OUTPUT = True
COMBINED_SMILES_FILE = r"D:\TNBC\bharati_dock\all_compounds.smi"


# ═════════════════════════════════════════════════════════════════════════════
#  METHOD 1 — via openbabel Python bindings (openbabel-wheel)
# ═════════════════════════════════════════════════════════════════════════════
def convert_with_pybel(input_folder, output_folder, combined):
    try:
        from openbabel import pybel
    except ImportError:
        return False, "openbabel Python bindings not found"

    os.makedirs(output_folder, exist_ok=True)
    files   = sorted(Path(input_folder).glob("*.pdbqt"))
    success = 0
    failed  = []
    combined_records = []

    print(f"\n  Using : OpenBabel Python bindings (pybel)")
    print(f"  Found : {len(files)} .pdbqt file(s)\n")

    for fpath in files:
        try:
            mols = list(pybel.readfile("pdbqt", str(fpath)))
            if not mols:
                failed.append((fpath.name, "no molecule parsed"))
                continue

            mol = mols[0]
            mol.removeh()                          # strip explicit H for cleaner SMILES
            smi = mol.write("smi").strip()
            # pybel.write returns "SMILES  name\n" — split to get pure SMILES
            smi_only = smi.split()[0] if smi else ""

            if not smi_only:
                failed.append((fpath.name, "empty SMILES generated"))
                continue

            stem = fpath.stem
            record = f"{smi_only}\t{stem}"
            combined_records.append(record)

            if not combined:
                out_path = Path(output_folder) / (stem + ".smi")
                out_path.write_text(record + "\n")
                print(f"  [OK]  {fpath.name}  →  {out_path.name}")
                print(f"        SMILES: {smi_only[:80]}{'...' if len(smi_only)>80 else ''}")
            else:
                print(f"  [OK]  {fpath.name}")
                print(f"        SMILES: {smi_only[:80]}{'...' if len(smi_only)>80 else ''}")

            success += 1

        except Exception as e:
            failed.append((fpath.name, str(e)))
            print(f"  [FAIL] {fpath.name} — {e}")

    # Write combined file
    if combined and combined_records:
        comb_path = Path(COMBINED_SMILES_FILE)
        comb_path.parent.mkdir(parents=True, exist_ok=True)
        comb_path.write_text("\n".join(combined_records) + "\n")
        print(f"\n  Combined SMI saved to: {COMBINED_SMILES_FILE}")

    return True, (success, failed, len(files))


# ═════════════════════════════════════════════════════════════════════════════
#  METHOD 2 — via obabel command-line tool
# ═════════════════════════════════════════════════════════════════════════════
def obabel_available():
    try:
        result = subprocess.run(
            ["obabel", "--version"],
            capture_output=True, text=True, timeout=10
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def convert_with_obabel_cli(input_folder, output_folder, combined):
    if not obabel_available():
        return False, "obabel command-line tool not found in PATH"

    os.makedirs(output_folder, exist_ok=True)
    files   = sorted(Path(input_folder).glob("*.pdbqt"))
    success = 0
    failed  = []
    combined_records = []

    print(f"\n  Using : obabel command-line tool")
    print(f"  Found : {len(files)} .pdbqt file(s)\n")

    for fpath in files:
        stem     = fpath.stem
        out_path = Path(output_folder) / (stem + ".smi")

        try:
            cmd = [
                "obabel",
                str(fpath),
                "-isdf" if fpath.suffix.lower() == ".sdf" else "-ipdbqt",
                "-osmi",
                "-O", str(out_path),
                "--gen2D",      # generate 2D coordinates (improves SMILES quality)
                "-h",           # add hydrogens then strip for canonical SMILES
            ]
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=60
            )

            if result.returncode != 0 or not out_path.exists():
                err = result.stderr.strip() or "unknown error"
                failed.append((fpath.name, err))
                print(f"  [FAIL] {fpath.name} — {err}")
                continue

            smi_content = out_path.read_text().strip()
            if not smi_content:
                failed.append((fpath.name, "empty output"))
                print(f"  [FAIL] {fpath.name} — empty SMILES output")
                continue

            smi_line = smi_content.split("\n")[0]
            smi_only = smi_line.split()[0] if smi_line else ""

            print(f"  [OK]  {fpath.name}  →  {out_path.name}")
            print(f"        SMILES: {smi_only[:80]}{'...' if len(smi_only)>80 else ''}")

            combined_records.append(f"{smi_only}\t{stem}")
            success += 1

        except subprocess.TimeoutExpired:
            failed.append((fpath.name, "timeout"))
            print(f"  [FAIL] {fpath.name} — conversion timed out")
        except Exception as e:
            failed.append((fpath.name, str(e)))
            print(f"  [FAIL] {fpath.name} — {e}")

    # Write combined file
    if combined and combined_records:
        comb_path = Path(COMBINED_SMILES_FILE)
        comb_path.parent.mkdir(parents=True, exist_ok=True)
        comb_path.write_text("\n".join(combined_records) + "\n")
        print(f"\n  Combined SMI saved to: {COMBINED_SMILES_FILE}")

    return True, (success, failed, len(files))


# ═════════════════════════════════════════════════════════════════════════════
#  REPORT
# ═════════════════════════════════════════════════════════════════════════════
def print_report(success, failed, total):
    print("\n" + "=" * 70)
    print(f"  CONVERSION SUMMARY")
    print(f"  Total files : {total}")
    print(f"  Succeeded   : {success}")
    print(f"  Failed      : {len(failed)}")
    if failed:
        print("\n  Failed files:")
        for fname, reason in failed:
            print(f"    ✗  {fname}  ({reason})")
    print("=" * 70)

    if COMBINED_OUTPUT:
        print(f"\n  All SMILES combined in : {COMBINED_SMILES_FILE}")
    else:
        print(f"\n  Individual .smi files  : {OUTPUT_FOLDER}")
    print()


# ═════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═════════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 70)
    print("  PDBQT → SMI Converter")
    print(f"  Input  : {INPUT_FOLDER}")
    if COMBINED_OUTPUT:
        print(f"  Output : {COMBINED_SMILES_FILE}  (combined)")
    else:
        print(f"  Output : {OUTPUT_FOLDER}  (individual files)")
    print("=" * 70)

    # Validate input folder
    if not os.path.isdir(INPUT_FOLDER):
        print(f"\nERROR: Input folder not found:\n  {INPUT_FOLDER}")
        sys.exit(1)

    files = list(Path(INPUT_FOLDER).glob("*.pdbqt"))
    if not files:
        print(f"\nNo .pdbqt files found in:\n  {INPUT_FOLDER}")
        sys.exit(0)

    # Try Method 1: pybel (openbabel Python bindings)
    ok, result = convert_with_pybel(INPUT_FOLDER, OUTPUT_FOLDER, COMBINED_OUTPUT)
    if ok:
        success, failed, total = result
        print_report(success, failed, total)
        return

    print(f"\n  [INFO] pybel not available ({result})")
    print("  [INFO] Trying obabel command-line tool...\n")

    # Try Method 2: obabel CLI
    ok, result = convert_with_obabel_cli(INPUT_FOLDER, OUTPUT_FOLDER, COMBINED_OUTPUT)
    if ok:
        success, failed, total = result
        print_report(success, failed, total)
        return

    print(f"\n  [INFO] obabel CLI not available ({result})")

    # Neither method available — guide the user
    print("\n" + "=" * 70)
    print("  OpenBabel NOT FOUND")
    print("=" * 70)
    print("""
  This script requires OpenBabel for PDBQT → SMILES conversion.
  Please install it using ONE of the following methods:

  ── Option A (Recommended — Python package) ──────────────────────────
      pip install openbabel-wheel

  ── Option B (Full installer for Windows) ────────────────────────────
      1. Download from: https://github.com/openbabel/openbabel/releases
      2. Install and add to PATH
      3. Restart terminal and re-run this script

  After installation, re-run:
      python pdbqt_to_smi.py
""")
    sys.exit(1)


if __name__ == "__main__":
    main()