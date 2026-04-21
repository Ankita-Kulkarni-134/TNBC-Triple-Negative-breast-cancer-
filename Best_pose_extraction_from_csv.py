import csv
import os

# Paths
INPUT_CSV  = r"D:\TNBC\bharati_dock\Result_600_extrcated\600MW_BG_greater than -8.5 - Copy.csv"
OUTPUT_CSV = r"D:\TNBC\bharati_dock\Result_600_extrcated\unique_best_VINA.csv"

KEY_COL     = "Filename"
COMPARE_COL = "VINA_Result_1"

def main():
    if not os.path.isfile(INPUT_CSV):
        print("ERROR: Input file not found:")
        print("  " + INPUT_CSV)
        return

    print("=" * 65)
    print("  Extract Unique Filenames - Best VINA_Result_1")
    print("  Input  : " + INPUT_CSV)
    print("  Output : " + OUTPUT_CSV)
    print("=" * 65)

    best_rows  = {}
    total_rows = 0
    skipped    = 0

    with open(INPUT_CSV, newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)

        if KEY_COL not in reader.fieldnames:
            print("ERROR: Column '" + KEY_COL + "' not found in CSV.")
            print("  Available columns: " + str(reader.fieldnames))
            return
        if COMPARE_COL not in reader.fieldnames:
            print("ERROR: Column '" + COMPARE_COL + "' not found in CSV.")
            print("  Available columns: " + str(reader.fieldnames))
            return

        fieldnames = reader.fieldnames

        for row in reader:
            total_rows += 1
            filename = row[KEY_COL].strip()

            try:
                score = float(row[COMPARE_COL])
            except (ValueError, TypeError):
                skipped += 1
                continue

            if filename not in best_rows:
                best_rows[filename] = (score, row)
            else:
                current_best_score, _ = best_rows[filename]
                if score > current_best_score:
                    best_rows[filename] = (score, row)

    unique_count = len(best_rows)

    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)

    with open(OUTPUT_CSV, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for fname in sorted(best_rows.keys()):
            score, row = best_rows[fname]
            writer.writerow(row)

    print("")
    print("  Total rows read        : " + str(total_rows))
    print("  Rows skipped (errors)  : " + str(skipped))
    print("  Unique filenames kept  : " + str(unique_count))
    print("  Duplicate rows removed : " + str(total_rows - skipped - unique_count))
    print("")
    print("  Output saved to:")
    print("  " + OUTPUT_CSV)
    print("=" * 65)

    print("")
    print("  Top 10 entries (best VINA_Result_1):")
    print("")
    print("  {:<30} {:>6}  {:>14}".format("Filename", "Model", "VINA_Result_1"))
    print("  " + "-" * 55)

    sorted_results = sorted(best_rows.items(), key=lambda x: x[1][0], reverse=True)

    for fname, (score, row) in sorted_results[:10]:
        model = row.get("Model", "?")
        print("  {:<30} {:>6}  {:>14.1f}".format(fname, model, score))

    print("")

if __name__ == "__main__":
    main()