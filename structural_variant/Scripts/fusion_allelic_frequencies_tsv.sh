#!/usr/bin/env bash

# ============================
# CONFIGURATION
# ============================

files=(
    "/beegfs/data/fblanchard/TE/structural_variant/Results/S10_Wb_noTE_nolow_SV/graphaligner/allelic_frequencies.tsv"
    "/beegfs/data/fblanchard/TE/structural_variant/Results/S31_Wb_noTE_nolow_SV/graphaligner/allelic_frequencies.tsv"
    "/beegfs/data/fblanchard/TE/structural_variant/Results/S73_Wb_noTE_nolow_SV/graphaligner/allelic_frequencies.tsv"
    "/beegfs/data/fblanchard/TE/structural_variant/Results/S100_Wb_noTE_nolow_SV/graphaligner/allelic_frequencies.tsv"
)

vcf_files=(
    "/beegfs/data/fblanchard/TE/structural_variant/Results/S10_Wb_noTE_nolow_SV/graphaligner/S10_Wb_noTE_nolow_SV_graphaligner.vcf"
    "/beegfs/data/fblanchard/TE/structural_variant/Results/S31_Wb_noTE_nolow_SV/graphaligner/S31_Wb_noTE_nolow_SV_graphaligner.vcf"
    "/beegfs/data/fblanchard/TE/structural_variant/Results/S73_Wb_noTE_nolow_SV/graphaligner/S73_Wb_noTE_nolow_SV_graphaligner.vcf"
    "/beegfs/data/fblanchard/TE/structural_variant/Results/S100_Wb_noTE_nolow_SV/graphaligner/S100_Wb_noTE_nolow_SV_graphaligner.vcf"
)

output="/beegfs/data/fblanchard/TE/structural_variant/Results/resume.tsv"


# ============================
# VERIFICATIONS
# ============================

if [ ${#files[@]} -eq 0 ] || [ ${#vcf_files[@]} -eq 0 ]; then
  echo "❌ TSV ou VCF manquants."
  exit 1
fi

if [ ${#files[@]} -ne ${#vcf_files[@]} ]; then
  echo "❌ Le nombre de TSV ne correspond pas au nombre de VCF."
  exit 1
fi


# ============================
# ETAPE 1 : EXTRACTION DES CHR:POS À EXCLURE
# ============================

echo "⏳ Filtrage des sites avec ≥3 N consécutifs dans REF/ALT..."

tmp_exclude_ids=$(mktemp)

for vcf in "${vcf_files[@]}"; do
    awk -v FS='\t' '
        $0 !~ /^#/ {
            id = $1 ":" $2
            if ($4 ~ /NNN/ || $5 ~ /NNN/) print id
        }
    ' "$vcf"
done | sort -u > "$tmp_exclude_ids"

echo "➡️  $(wc -l < "$tmp_exclude_ids") sites exclus."


# ============================
# ETAPE 2 : EN-TÊTE DU FICHIER FINAL
# ============================

header="CHR_POS"
for f in "${files[@]}"; do
    name=$(basename "$f" .tsv)
    header+="\t${name}"
done
header+="\tTYPE"

echo -e "$header" > "$output"


# ============================
# ETAPE 3 : FUSION + FILTRES + TYPE INS/DEL
# ============================

echo "⏳ Fusion des TSV + classification INS/DEL..."

awk -v FS='\t' -v OFS='\t' \
    -v exclude="$tmp_exclude_ids" \
    -v vcf_list="${vcf_files[*]}" '
BEGIN {
    # Charger CHR:POS exclus
    while ((getline line < exclude) > 0)
        excluded[line] = 1

    # Charger type INS/DEL depuis les VCF
    split(vcf_list, vcf_arr, " ")

    for (v in vcf_arr) {
        file = vcf_arr[v]
        while ((getline < file) > 0) {
            if ($0 ~ /^#/) continue

            id = $1 ":" $2
            ref = $4
            alt = $5

            if (type[id] != "") continue

            if (ref ~ /^[ACGTN]$/ && length(alt) > 1)
                type[id] = "INS"
            else if (alt ~ /^[ACGTN]$/ && length(ref) > 1)
                type[id] = "DEL"
            else
                type[id] = "NA"
        }
        close(file)
    }
}
{
    # Utiliser colonnes 1 et 2 du TSV comme clé
    id = $1 ":" $2

    # Si exclu → ignorer
    if (excluded[id] == 1) next

    val = $4
    data[FILENAME][id] = val
    all_ids[id] = 1
}
END {
    n = ARGC - 1
    PROCINFO["sorted_in"] = "@ind_str_asc"

    for (id in all_ids) {

        zero_count = 0
        line = id

        for (i = 1; i <= n; i++) {
            val = data[ARGV[i]][id]
            if (val == "") val = 0
            if (val == 0) zero_count++
            line = line "\t" val
        }

        # Si toutes les lignes = 0 → ignorer
        if (zero_count == n) continue

        # Ajouter le type INS / DEL / NA
        if (type[id] == "") t = "NA"
        else t = type[id]

        print line "\t" t
    }
}' "${files[@]}" >> "$output"

echo "✅ Fichier final créé : $output"
