#!/usr/bin/env python3

import gzip
import sys

# Fichiers VCF
vcf_files = [
    "/beegfs/data/fblanchard/TE/structural_variant/Results/S10_Wb_noTE_nolow_SV/graphaligner/S10_Wb_noTE_nolow_SV_graphaligner.vcf",
    "/beegfs/data/fblanchard/TE/structural_variant/Results/S31_Wb_noTE_nolow_SV/graphaligner/S31_Wb_noTE_nolow_SV_graphaligner.vcf",
    "/beegfs/data/fblanchard/TE/structural_variant/Results/S73_Wb_noTE_nolow_SV/graphaligner/S73_Wb_noTE_nolow_SV_graphaligner.vcf",
    "/beegfs/data/fblanchard/TE/structural_variant/Results/S100_Wb_noTE_nolow_SV/graphaligner/S100_Wb_noTE_nolow_SV_graphaligner.vcf"
]

table_file = "/beegfs/data/fblanchard/TE/structural_variant/Results/positions_proches_clusters2plus.tsv"   # ton tableau d'entrée
output_fasta = "/beegfs/data/fblanchard/TE/structural_variant/Results/SV_sequences.fasta"

# ---------------------------
# INDEXAGE DES VCF PAR CHROM+POS
# ---------------------------

index = {}  # (chrom, pos) → {"REF":..., "ALT":...}

def open_vcf(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    return open(filename, "r")

for vcf in vcf_files:
    with open_vcf(vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")

            chrom = cols[0]
            pos   = int(cols[1])
            ref   = cols[3]
            alt   = cols[4]

            index[(chrom, pos)] = {"REF": ref, "ALT": alt}

print(f"Index VCF : {len(index)} variants chargés.")


# ---------------------------
# PARCOURS DU TABLEAU
# ---------------------------

with open(table_file) as t, open(output_fasta, "w") as out:
    for line in t:
        if line.strip() == "":
            continue

        cols = line.strip().split("\t")
        chrom = cols[0]
        pos = int(cols[1])
        svtype = cols[6]  # DEL ou INS

        key = (chrom, pos)

        if key not in index:
            print(f"[WARNING] Variante {chrom}:{pos} introuvable dans les VCF.")
            continue

        ref = index[key]["REF"]
        alt = index[key]["ALT"]

        # DEL → REF | INS → ALT
        seq = ref if svtype == "DEL" else alt

        header = f">{chrom}_{pos}"
        out.write(header + "\n")
        out.write(seq + "\n")

print(f"FASTA écrit : {output_fasta}")
