#!/usr/bin/env python3
"""
plasmid_localize_amr_genes.py

Joins Abricate CARD + PlasmidFinder hits with MOB-suite's contig-to-plasmid
assignments, to see whether each hit sits on a chromosome-associated contig or
a MOB-reconstructed plasmid unit.

The CARD gene -> resistance class mapping below is a starting point, not a
validated mapping. The script writes out every CARD gene it couldn't map
(_unmapped_genes.csv); check that, and spot-check the mapped genes against
Kleborate's calls in character_matrix.csv, before using resistance_class for
anything.

genome_id is the filename minus .tab (files are <GCF_id>.tab, same as the
genome_id column in character_matrix.csv).

Usage:
    python plasmid_localize_amr_genes.py \\
        --abricate-dir abricate_out \\
        --mob-contig-table contig_to_cluster.tsv \\
        --out-prefix phase10_out/amr_plasmid_localization \\
        [--hits-are-molecule]

Inputs:
    --abricate-dir/card/<GCF_id>.tab           Abricate CARD output
    --abricate-dir/plasmidfinder/<GCF_id>.tab  Abricate PlasmidFinder output
    --mob-contig-table                         MOB-suite contig_to_cluster.tsv
                                               (genome_id, molecule_type, primary_cluster_id,
                                               secondary_cluster_id, contig_id, ...)

Outputs:
    <out-prefix>.csv                     one row per Abricate hit, with molecule_type and
                                         primary_cluster_id joined in where matched
    <out-prefix>_unmatched_contigs.csv   hits whose contig_id isn't in the MOB table at all
    <out-prefix>_unmapped_genes.csv      CARD genes that didn't map to any resistance class
    <out-prefix>_qc.txt                  summary counts, read this before trusting the rest
"""

import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

# CARD gene name -> resistance class. First matching pattern wins, so order
# matters. Matched case-insensitively against Abricate's GENE column.
CARD_TO_CLASS_PATTERNS = [
    # carbapenemases (this CARD version uses bare allele names like "KPC-2", no "bla" prefix)
    (r"^KPC-", "Bla_Carb_acquired"),
    (r"^NDM-", "Bla_Carb_acquired"),
    (r"^IMP-", "Bla_Carb_acquired"),
    (r"^VIM-", "Bla_Carb_acquired"),
    # only the OXA-48-like alleles are carbapenemases; other OXAs (OXA-1, -2, -9, -10...)
    # aren't. OXA-9 is penicillin-only per CARD's RESISTANCE field, it just co-occurs
    # with KPC on shared plasmids
    (r"^OXA-(48|181|232|244|245|204|162|163|517|519|544|796|827|896|926|933|"
      r"1136|1167|1175|1201|1306)$", "Bla_Carb_acquired"),
    # ESBLs
    (r"^CTX-M-", "Bla_ESBL_acquired"),
    (r"^VEB-", "Bla_ESBL_acquired"),
    (r"^PER-", "Bla_ESBL_acquired"),
    # GES is split like OXA (most are ESBL, GES-5 and a few hydrolyse carbapenems),
    # so it's left unclassified until alleles are checked
    (r"^GES-", "GES_family_needs_allele_check"),
    # remaining OXAs: narrow-spectrum
    (r"^OXA-", "Bla_acquired"),
    # acquired AmpC has no slot in the locked 3-way Bla_* split; kept as its own
    # label instead of reopening a frozen Phase 2 decision
    (r"^DHA-|^CMY-|^FOX-|^MOX-|^ACC-|^LAP-", "AmpC_acquired_not_in_locked_set"),
    # LEN/OKP/KLUC are intrinsic chromosomal beta-lactamases (like SHV-1), not acquired
    (r"^LEN-|^OKP-|^KLUC-", "intrinsic_chromosomal_not_acquired"),
    # SHV and TEM: narrow vs ESBL vs BLI-resistant depends on specific amino-acid
    # substitutions, not allele names (Liu et al. 2024). Kleborate v2.4.1+ already does
    # this properly in character_matrix.csv, so look the class up there instead
    (r"^SHV-", "SHV_family_see_Kleborate_call"),
    (r"^TEM-", "TEM_family_see_Kleborate_call"),
    # intrinsic K. pneumoniae efflux/porin/regulator/LPS genes: present in ~all genomes,
    # not acquired resistance. Exact names only so they can't swallow other genes
    (r"^(OmpA|CRP|H-NS|msbA|eptA|eptB|LptD|cpxA|MdtQ|ArnT|marA|ramA|emrR|"
     r"ugd|basS|LpsB|ArmR|rsmA|CrcB|abeM|MexD|AcrF|bacA|kdpE|gadW|gadX|evgA|evgS|"
     r"acrB|acrD|emrA|emrB|emrK|emrY|mdt[A-H]|"
     r"Klebsiella_pneumoniae_(acrA|KpnE|KpnF|KpnG|KpnH|OmpK37)|"
     r"Escherichia_coli_(acrA|emrE))$", "intrinsic_chromosomal_efflux_regulatory"),
    # class 1 integron / quaternary ammonium markers
    (r"^qacE(delta1)?$|^qacL$", "disinfectant_qac_not_in_locked_set"),
    # bleomycin resistance (BRP(MBL) is the NDM-linked marker, not an actual MBL)
    (r"^BRP\(MBL\)$|^BLMT$", "bleomycin_not_in_locked_set"),
    # streptothricin acetyltransferases, not an aminoglycoside class in the locked set
    (r"^SAT-\d", "streptothricin_not_in_locked_set"),
    (r"^Mrx$", "MLS_acquired"),                 # (verify) mph/mrx macrolide locus
    (r"^EstT$", "Rif_acquired"),                # rifamycin-inactivating esterase
    (r"^cfrA$", "Phe_acquired"),                # (verify) Cfr: phenicol + others, 3 hits
    # SFO-1 / FRI-1 / RSA1-1: CARD's drug class conflicts with the literature, so left out
    # of the locked Bla_* classes until checked per allele (28 hits)
    (r"^(SFO-1|FRI-1|RSA1-1)$", "Bla_family_needs_allele_check"),
    (r"^(SCO-1|CARB-3|CARB-58)$", "Bla_acquired"),  # (verify) narrow-spectrum
    (r"^(SED-2|EC-18)$", "chromosomal_or_other_bla_not_in_locked_set"),  # (verify)
    # oqxAB is near-universal in K. pneumoniae and would inflate Flq, so it's split out
    (r"^oqx[AB]$", "oqxAB_efflux_not_in_locked_set"),
    # FosA6 is the chromosomal K. pneumoniae fosA (~88% of genomes, 7,356 hits), not
    # acquired resistance. Other alleles (FosA3/4/5/7/8) stay Fcyn via the generic fosA rule
    (r"^fosA[567]$", "fosA_chromosomal_not_in_locked_set"),
    (r"aac\(6.\)-Ib-cr", "Flq_acquired"),
    (r"aac\(6.?\)|aac\(3\)|aac\(2.?\)|aph\(|aphA\d|ant\(|aad[ab]|^rmt[A-Z]$|^armA$", "AGly_acquired"),
    (r"^mcr-", "Col_acquired"),
    (r"fosA", "Fcyn_acquired"),
    (r"qnr|oqx|aac\(6.\)-Ib-cr|gyrA|parC|^QepA", "Flq_acquired"),
    (r"erm|mph|msr|lnu|^Ere[AB]|^mef\(B\)|^linG", "MLS_acquired"),
    (r"catA|cml|floR|catB|^catII|^catIII|^cmx", "Phe_acquired"),
    (r"^arr-", "Rif_acquired"),
    (r"^sul[123]", "Sul_acquired"),
    (r"^tmexCD|^tet\(X\)|^tigecycline", "Tgc_acquired"),
    (r"^tet\(", "Tet_acquired"),
    (r"^tet[A-Z]$", "Tet_acquired"),
    (r"^dfrA|^dfrB|^dfrG|^dfr\d|DfrA\d", "Tmt_acquired"),
]
COMPILED = [(re.compile(p, re.IGNORECASE), c) for p, c in CARD_TO_CLASS_PATTERNS]


def map_gene_to_class(gene_name: str):
    for pattern, cls in COMPILED:
        if pattern.search(gene_name):
            return cls
    return None


def genome_id_from_path(path: str) -> str:
    # <GCF_id>.tab -> <GCF_id>
    base = os.path.basename(path)
    if base.endswith(".tab"):
        return base[: -len(".tab")]
    raise ValueError(f"Expected a .tab file, got: {path}")


def read_abricate_tab(path: str, database_label: str):
    # one dict per hit row; a header-only file (no hits) yields nothing, which is fine
    genome_id = genome_id_from_path(path)
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            return
        for row in reader:
            yield {
                "genome_id": genome_id,
                "database": database_label,
                "contig_id": row.get("SEQUENCE", ""),
                "gene": row.get("GENE", ""),
                "pct_coverage": row.get("%COVERAGE", ""),
                "pct_identity": row.get("%IDENTITY", ""),
                "product": row.get("PRODUCT", ""),
                "resistance_field": row.get("RESISTANCE", ""),
            }


def _first_token(s: str) -> str:
    # MOB's contig_id holds the full FASTA header ('NZ_CP006659.2 Klebsiella ...');
    # Abricate's SEQUENCE only keeps the first word. Compare on the first token.
    return s.split()[0] if s.split() else s


def load_mob_contig_table(path: str):
    # (genome_id, first token of contig_id) -> molecule_type + cluster ids.
    # The table has both genome_id (clean accession) and sample_id (raw input
    # filename); genome_id is the one that matches Abricate's.
    expected_cols = {"genome_id", "molecule_type", "primary_cluster_id",
                     "secondary_cluster_id", "contig_id"}
    lookup = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = expected_cols - set(reader.fieldnames or [])
        if missing:
            sys.exit(
                f"ERROR: MOB contig table is missing expected columns: {missing}\n"
                f"Actual columns found: {reader.fieldnames}\n"
                f"Refusing to proceed with a silent column mismatch."
            )
        for row in reader:
            key = (row["genome_id"], _first_token(row["contig_id"]))
            lookup[key] = {
                "molecule_type": row["molecule_type"],
                "primary_cluster_id": row["primary_cluster_id"],
                "secondary_cluster_id": row["secondary_cluster_id"],
            }
    return lookup


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--abricate-dir", required=True,
                        help="Directory containing card/ and plasmidfinder/ subfolders of .tab files")
    parser.add_argument("--mob-contig-table", required=True,
                        help="MOB-suite's contig_to_cluster.tsv")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument(
        "--hits-are-molecule", action="store_true",
        help=(
            "Skip the MOB contig-table join and infer molecule type from the contig_id "
            "string itself ('plasmid' or 'chromosome' in the name, else 'unknown'). "
            "Default: join every hit against --mob-contig-table by (genome_id, contig_id)."
        ),
    )
    args = parser.parse_args()

    card_files = sorted(glob.glob(os.path.join(args.abricate_dir, "card", "*.tab")))
    pf_files = sorted(glob.glob(os.path.join(args.abricate_dir, "plasmidfinder", "*.tab")))
    if not card_files or not pf_files:
        sys.exit(
            f"ERROR: found {len(card_files)} card/*.tab and {len(pf_files)} "
            f"plasmidfinder/*.tab under {args.abricate_dir}. Expected 8318 each. "
            f"Check --abricate-dir before proceeding."
        )

    mob_lookup = None
    if not args.hits_are_molecule:
        mob_lookup = load_mob_contig_table(args.mob_contig_table)

    all_hits = []
    for path in card_files:
        all_hits.extend(read_abricate_tab(path, "card"))
    for path in pf_files:
        all_hits.extend(read_abricate_tab(path, "plasmidfinder"))

    unmatched_contigs = []
    unmapped_genes = Counter()
    molecule_counts = Counter()

    for hit in all_hits:
        # molecule type
        if args.hits_are_molecule:
            cid_lower = hit["contig_id"].lower()
            if "plasmid" in cid_lower:
                hit["molecule_type"] = "plasmid"
            elif "chromosome" in cid_lower:
                hit["molecule_type"] = "chromosome"
            else:
                hit["molecule_type"] = "unknown"
            hit["primary_cluster_id"] = ""
            hit["secondary_cluster_id"] = ""
        else:
            key = (hit["genome_id"], _first_token(hit["contig_id"]))
            mob_row = mob_lookup.get(key)
            if mob_row is None:
                hit["molecule_type"] = "UNMATCHED"
                hit["primary_cluster_id"] = ""
                hit["secondary_cluster_id"] = ""
                unmatched_contigs.append(dict(hit))
            else:
                hit["molecule_type"] = mob_row["molecule_type"]
                hit["primary_cluster_id"] = mob_row["primary_cluster_id"]
                hit["secondary_cluster_id"] = mob_row["secondary_cluster_id"]
        molecule_counts[hit["molecule_type"]] += 1

        # resistance class (CARD hits only)
        if hit["database"] == "card":
            cls = map_gene_to_class(hit["gene"])
            hit["resistance_class"] = cls if cls else ""
            if cls is None:
                unmapped_genes[hit["gene"]] += 1
        else:
            hit["resistance_class"] = ""

    # main output
    out_cols = ["genome_id", "database", "gene", "resistance_class",
                "contig_id", "molecule_type", "primary_cluster_id",
                "secondary_cluster_id", "pct_coverage", "pct_identity",
                "product", "resistance_field"]
    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    with open(f"{args.out_prefix}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols)
        w.writeheader()
        for hit in all_hits:
            w.writerow({k: hit.get(k, "") for k in out_cols})

    with open(f"{args.out_prefix}_unmatched_contigs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols)
        w.writeheader()
        for hit in unmatched_contigs:
            w.writerow({k: hit.get(k, "") for k in out_cols})

    with open(f"{args.out_prefix}_unmapped_genes.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "hit_count"])
        for gene, count in unmapped_genes.most_common():
            w.writerow([gene, count])

    n_genomes = len({h["genome_id"] for h in all_hits})
    with open(f"{args.out_prefix}_qc.txt", "w") as f:
        f.write(f"Total hit rows: {len(all_hits)}\n")
        f.write(f"Distinct genomes with >=1 hit (either DB): {n_genomes}\n")
        f.write(f"CARD files read: {len(card_files)}\n")
        f.write(f"PlasmidFinder files read: {len(pf_files)}\n")
        f.write(f"Mode: {'hits-are-molecule (direct)' if args.hits_are_molecule else 'MOB contig-table join'}\n")
        f.write(f"Molecule-type breakdown: {dict(molecule_counts)}\n")
        if not args.hits_are_molecule:
            f.write(f"UNMATCHED hits (contig_id not found in MOB table): {len(unmatched_contigs)}\n")
            if unmatched_contigs:
                f.write(
                    "  >>> Non-zero UNMATCHED count. This means some Abricate contig_id "
                    "values don't appear in MOB's contig_to_cluster.tsv at all -- check "
                    "whether contig-ID naming differs between Abricate's SEQUENCE column "
                    "and MOB's contig_id column before trusting molecule_type broadly.\n"
                )
        f.write(f"CARD genes with NO resistance-class mapping: {len(unmapped_genes)} distinct genes, "
                f"{sum(unmapped_genes.values())} total hit rows\n")
        if unmapped_genes:
            f.write("  >>> Check _unmapped_genes.csv before trusting resistance_class downstream.\n")
            f.write("  Top unmapped genes: " +
                    ", ".join(f"{g}({c})" for g, c in unmapped_genes.most_common(10)) + "\n")

    print(f"Done. Wrote {args.out_prefix}.csv ({len(all_hits)} rows). "
          f"Read the _qc.txt file before using this output for anything.")


if __name__ == "__main__":
    main()