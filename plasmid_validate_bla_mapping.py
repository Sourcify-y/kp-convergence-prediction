#!/usr/bin/env python3
"""
plasmid_validate_bla_mapping.py

For each SHV/TEM/OXA/AmpC-family/LEN/OKP CARD allele, check how often the
genome also has Kleborate's Bla_acquired / Bla_ESBL_acquired / Bla_Carb_acquired
call set to 1 in character_matrix.csv. This gives an empirical answer to which
class each allele belongs in, instead of trusting allele-name regex guesses.

Usage:
    python plasmid_validate_bla_mapping.py \\
        --localization-csv test_out/localization_test.csv \\
        --character-matrix phase3_out/character_matrix.csv \\
        --out concordance_report.csv

Only CARD hits whose gene starts with SHV-, TEM-, OXA-, DHA-, CMY-, FOX-,
MOX-, ACC-, LEN- or OKP- are used (the groups flagged as unvalidated in
plasmid_localize_amr_genes.py). A clean allele agrees ~100% with exactly one
Bla_* column; mixed agreement means either the allele really is variable in
effect or Kleborate and CARD disagree for some other reason worth a look.
"""

import argparse
import csv
import re
from collections import defaultdict

FLAGGED_RE = re.compile(r"^(SHV|TEM|OXA|DHA|CMY|FOX|MOX|ACC|LEN|OKP)-")
BLA_CLASSES = ("Bla_acquired", "Bla_ESBL_acquired", "Bla_Carb_acquired")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--localization-csv", required=True)
    parser.add_argument("--character-matrix", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    # genome_id -> {Bla class column: "0"/"1"}
    kleb_calls = {}
    with open(args.character_matrix, newline="") as f:
        reader = csv.DictReader(f)
        needed = {"genome_id", "Bla_acquired", "Bla_ESBL_acquired", "Bla_Carb_acquired"}
        missing = needed - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"character_matrix.csv missing columns: {missing}\n"
                             f"Actual columns: {reader.fieldnames}")
        for row in reader:
            kleb_calls[row["genome_id"]] = {
                "Bla_acquired": row["Bla_acquired"],
                "Bla_ESBL_acquired": row["Bla_ESBL_acquired"],
                "Bla_Carb_acquired": row["Bla_Carb_acquired"],
            }

    # allele -> {Bla class: number of its hits where Kleborate's call was 1}
    cooc = defaultdict(lambda: defaultdict(int))
    allele_n = defaultdict(int)

    with open(args.localization_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["database"] != "card":
                continue
            gene = row["gene"]
            if not FLAGGED_RE.match(gene):
                continue
            genome_id = row["genome_id"]
            calls = kleb_calls.get(genome_id)
            if calls is None:
                continue
            allele_n[gene] += 1
            for cls in BLA_CLASSES:
                if calls[cls] == "1":
                    cooc[gene][cls] += 1

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["allele", "n_genomes", "pct_Bla_acquired",
                    "pct_Bla_ESBL_acquired", "pct_Bla_Carb_acquired",
                    "suggested_class"])
        for gene in sorted(allele_n, key=lambda g: -allele_n[g]):
            n = allele_n[gene]
            pcts = {}
            for cls in BLA_CLASSES:
                pcts[cls] = cooc[gene][cls] / n
            best_cls = max(pcts, key=pcts.get)
            best_pct = pcts[best_cls]
            suggestion = best_cls if best_pct >= 0.9 else f"UNCLEAR (best={best_cls} @ {best_pct:.0%})"
            w.writerow([gene, n, f"{pcts['Bla_acquired']:.2%}",
                        f"{pcts['Bla_ESBL_acquired']:.2%}",
                        f"{pcts['Bla_Carb_acquired']:.2%}", suggestion])

    print(f"Wrote {args.out}. Read it top to bottom -- alleles with a clean "
          f">=90% match to one class are safe to hardcode; UNCLEAR ones need "
          f"a manual look (or may indicate Kleborate/CARD genuinely disagree "
          f"for that allele, which is itself worth knowing).")


if __name__ == "__main__":
    main()