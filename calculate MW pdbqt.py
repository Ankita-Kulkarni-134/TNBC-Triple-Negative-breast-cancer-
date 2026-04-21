#!/usr/bin/env python3
"""
PDBQT Molecular Weight Calculator
===================================
Reads all .pdbqt files from:   D:\TNBC\bharati_dock\Result
Saves CSV results to:          D:\TNBC\bharati_dock\mw_results.csv

Run:
    python calculate_mw_pdbqt.py
"""

import os
import sys
from collections import defaultdict

# ── Paths ─────────────────────────────────────────────────────────────────────
INPUT_FOLDER  = r"D:\TNBC\bharati_dock\Result"
OUTPUT_CSV    = r"D:\TNBC\bharati_dock\mw_results.csv"

# ── Atomic masses (average) in g/mol ─────────────────────────────────────────
ATOMIC_MASSES = {
    'H':   1.008,  'HE':  4.003,  'LI':  6.941,  'BE':  9.012,
    'B':  10.811,  'C':  12.011,  'N':  14.007,  'O':  15.999,
    'F':  18.998,  'NE': 20.180,  'NA': 22.990,  'MG': 24.305,
    'AL': 26.982,  'SI': 28.086,  'P':  30.974,  'S':  32.065,
    'CL': 35.453,  'AR': 39.948,  'K':  39.098,  'CA': 40.078,
    'SC': 44.956,  'TI': 47.867,  'V':  50.942,  'CR': 51.996,
    'MN': 54.938,  'FE': 55.845,  'CO': 58.933,  'NI': 58.693,
    'CU': 63.546,  'ZN': 65.380,  'GA': 69.723,  'GE': 72.630,
    'AS': 74.922,  'SE': 78.971,  'BR': 79.904,  'KR': 83.798,
    'RB': 85.468,  'SR': 87.620,  'Y':  88.906,  'ZR': 91.224,
    'NB': 92.906,  'MO': 95.960,  'TC': 98.000,  'RU':101.070,
    'RH':102.906,  'PD':106.420,  'AG':107.868,  'CD':112.411,
    'IN':114.818,  'SN':118.710,  'SB':121.760,  'TE':127.600,
    'I': 126.904,  'XE':131.293,  'CS':132.905,  'BA':137.327,
    'LA':138.905,  'CE':140.116,  'PR':140.908,  'ND':144.242,
    'SM':150.360,  'EU':151.964,  'GD':157.250,  'TB':158.925,
    'DY':162.500,  'HO':164.930,  'ER':167.259,  'TM':168.934,
    'YB':173.054,  'LU':174.967,  'HF':178.490,  'TA':180.948,
    'W': 183.840,  'RE':186.207,  'OS':190.230,  'IR':192.217,
    'PT':195.084,  'AU':196.967,  'HG':200.592,  'TL':204.383,
    'PB':207.200,  'BI':208.980,  'U': 238.029,
}

# ── AutoDock PDBQT atom-type → element ───────────────────────────────────────
PDBQT_TYPE_TO_ELEMENT = {
    'C':  'C',   'A':  'C',
    'N':  'N',   'NA': 'N',   'NS': 'N',
    'O':  'O',   'OA': 'O',   'OS': 'O',
    'S':  'S',   'SA': 'S',
    'H':  'H',   'HD': 'H',   'HS': 'H',
    'P':  'P',   'F':  'F',   'I':  'I',
    'CL': 'CL',  'BR': 'BR',
    'FE': 'FE',  'ZN': 'ZN',  'MG': 'MG',
    'CA': 'CA',  'MN': 'MN',  'CU': 'CU',
    'NI': 'NI',  'CO': 'CO',  'K':  'K',
    'W':  'O',
    'G0': 'C',   'G1': 'C',   'G2': 'C',   'G3': 'C',  'G': 'C',
}


def element_from_pdbqt_type(ad_type: str):
    key = ad_type.strip().upper()
    mapped = PDBQT_TYPE_TO_ELEMENT.get(key)
    if mapped:
        return mapped.upper()
    for length in (2, 1):
        candidate = key[:length]
        if candidate in ATOMIC_MASSES:
            return candidate
    return None


def element_from_atom_name(atom_name: str):
    name = atom_name.strip().lstrip('0123456789')
    for length in (2, 1):
        candidate = name[:length].upper()
        if candidate in ATOMIC_MASSES:
            return candidate
    return None


def parse_pdbqt(filepath: str):
    total_mw    = 0.0
    elem_counts = defaultdict(int)
    warnings    = []
    atom_count  = 0

    with open(filepath, 'r') as fh:
        for lineno, line in enumerate(fh, 1):
            record = line[:6].strip().upper()
            if record not in ('ATOM', 'HETATM'):
                continue

            atom_count += 1
            element = None

            # 1) AutoDock type column (cols 77-79)
            if len(line) >= 78:
                ad_type = line[77:79].strip()
                if ad_type:
                    element = element_from_pdbqt_type(ad_type)

            # 2) Atom-name field (cols 13-16)
            if element is None and len(line) >= 16:
                element = element_from_atom_name(line[12:16])

            # 3) Standard PDB element column (cols 77-78)
            if element is None and len(line) >= 78:
                ef = line[76:78].strip().upper()
                if ef in ATOMIC_MASSES:
                    element = ef

            if element is None:
                raw = line[12:16].strip() if len(line) >= 16 else '?'
                warnings.append(f"  Line {lineno}: unknown element for atom '{raw}' — skipped")
                continue

            mass = ATOMIC_MASSES.get(element)
            if mass is None:
                warnings.append(f"  Line {lineno}: no mass for element '{element}' — skipped")
                continue

            total_mw += mass
            elem_counts[element] += 1

    return total_mw, dict(elem_counts), warnings, atom_count


def format_formula(element_counts: dict) -> str:
    order = ['C','H','N','O','S','P','F','CL','BR','I',
             'FE','ZN','MG','CA','MN','CU','NI','CO','NA','K']
    parts, seen = [], set()
    for el in order:
        if el in element_counts:
            parts.append(f"{el}:{element_counts[el]}")
            seen.add(el)
    for el in sorted(element_counts):
        if el not in seen:
            parts.append(f"{el}:{element_counts[el]}")
    return '  '.join(parts)


def main():
    if not os.path.isdir(INPUT_FOLDER):
        print(f"ERROR: Input folder not found:\n  {INPUT_FOLDER}")
        sys.exit(1)

    files = sorted(f for f in os.listdir(INPUT_FOLDER) if f.lower().endswith('.pdbqt'))

    if not files:
        print(f"No .pdbqt files found in:\n  {INPUT_FOLDER}")
        sys.exit(0)

    print("=" * 70)
    print("  PDBQT Molecular Weight Calculator")
    print(f"  Input  : {INPUT_FOLDER}")
    print(f"  Output : {OUTPUT_CSV}")
    print(f"  Files  : {len(files)} .pdbqt file(s) found")
    print("=" * 70)

    results = []

    for fname in files:
        fpath = os.path.join(INPUT_FOLDER, fname)
        mw, elem_counts, warns, n_atoms = parse_pdbqt(fpath)
        results.append((fname, mw, elem_counts, warns, n_atoms))

        print(f"\n  File     : {fname}")
        print(f"  Atoms    : {n_atoms}")
        print(f"  MW       : {mw:.3f} g/mol")
        if elem_counts:
            print(f"  Elements : {format_formula(elem_counts)}")
        if warns:
            print(f"  Warnings : {len(warns)}")
            for w in warns:
                print(w)

    # ── Summary table ─────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print(f"  {'Filename':<42} {'Atoms':>6}  {'MW (g/mol)':>12}")
    print("-" * 70)
    for fname, mw, _, _, n_atoms in results:
        print(f"  {fname:<42} {n_atoms:>6}  {mw:>12.3f}")
    print("=" * 70)

    # ── Write CSV ─────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    with open(OUTPUT_CSV, 'w') as csv_f:
        csv_f.write("Sr_No,Filename,Num_Atoms,Molecular_Weight_g_mol\n")
        for i, (fname, mw, _, _, n_atoms) in enumerate(results, 1):
            csv_f.write(f"{i},{fname},{n_atoms},{mw:.3f}\n")

    print(f"\n  CSV saved successfully to:\n  {OUTPUT_CSV}\n")


if __name__ == "__main__":
    main()