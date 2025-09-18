# paths
DB="processed/anml_aligned_recovered_dereplicated_filtered.txt"
PEST="priority_pests.txt"

# 1) extract unique binomial species from column 2 of the crabs db
#    - trims whitespace
#    - reduces any trinomial to genus + species
#    - drops entries like "Genus sp.", "Genus cf.", "Genus aff."
LC_ALL=C awk -F'\t' '
  {
    s=$2
    gsub(/^[[:space:]]+|[[:space:]]+$/,"",s)
    n=split(s,a," ")
    if (n>=2) {
      # skip non-binomial qualifiers
      if (tolower(a[2]) ~ /^(sp|cf|aff|nr)\.?$/) next
      printf("%s %s\n", tolower(a[1]), tolower(a[2]))
    }
  }
' "$DB" | sort -u > db_species.txt

# 2) normalize pest list to lowercased binomials (genus + species)
sed 's/^[[:space:]]\+//; s/[[:space:]]\+$//' "$PEST" \
  | sed 's/[[:space:]]\{1,\}/ /g' \
  | awk '
      { n=split(tolower($0),a," ");
        if (n>=2) printf("%s %s\n", a[1], a[2]);
      }
    ' \
  | sort -u > pest_clean.txt

# 3) compare sets
comm -12 pest_clean.txt db_species.txt > pests_present.txt
comm -23 pest_clean.txt db_species.txt > pests_missing.txt

# 4) quick summary
echo "present: $(wc -l < pests_present.txt)"
echo "missing: $(wc -l < pests_missing.txt)"

# 5) optional: counts per species remaining in db (based on column 2)
LC_ALL=C awk -F'\t' '
  {
    s=$2
    gsub(/^[[:space:]]+|[[:space:]]+$/,"",s)
    n=split(s,a," ")
    if (n>=2) {
      if (tolower(a[2]) ~ /^(sp|cf|aff|nr)\.?$/) next
      printf("%s %s\n", tolower(a[1]), tolower(a[2]))
    }
  }
' "$DB" \
| sort | uniq -c | sort -nr > db_species_counts.txt

# present: 21
# missing: 2
# the two missing:
# diaspidiotus perniciosus; typhlocyba pomaria

# 6) optional: counts only for pests that are present
grep -Fxf pests_present.txt db_species_counts.txt > pest_species_counts.txt || true




##############################################################################
#*                 make sure the missing pests are actually missing
##############################################################################
# paths
DB="processed/anml_aligned_recovered_OLD.txt"
PRE="processed/anml_aligned_recovered_dereplicated.txt"   # pre-filter
RAW="processed/anml_aligned_recovered.txt"                # pairwise-aligned output, pre-derep

# 0) minor fix: when reading a tabbed targets file with spaces in names, force IFS to tab
#    (your headings printed '== diaspidiotus ==' because read split on space)
# while IFS=$'\t' read -r canon pats; do ... done < targets.tsv

# 1) definitive present/absent in the final DB (exact, case-insensitive) for both species + common synonyms
for name in \
  "Diaspidiotus perniciosus" "Quadraspidiotus perniciosus" "Comstockaspis perniciosa" "Aspidiotus perniciosus" \
  "Typhlocyba pomaria" "Zonocyba pomaria" "Empoa pomaria"
do
  echo ">> final-db: $name"
  awk -F'\t' -v IGNORECASE=1 -v q="$name" '
    tolower($2)==tolower(q) || (NF>=10 && tolower($10)==tolower(q)) {print NR "\t" $0}
  ' "$DB" | head || true
done

# 2) if dias* is absent above, check whether it existed pre-filter or pre-derep
for F in "$PRE" "$RAW"; do
  echo "=== scan: $F"
  for name in "Diaspidiotus perniciosus" "Quadraspidiotus perniciosus" "Comstockaspis perniciosa" "Aspidiotus perniciosus"; do
    echo ">> $name"
    awk -F'\t' -v IGNORECASE=1 -v q="$name" '
      tolower($2)==tolower(q) || (NF>=10 && tolower($10)==tolower(q)) {print NR "\t" $0}
    ' "$F" | head || true
  done
done

# 3) broad sniff test for any 'pernicios*' token (helps catch unexpected genus assignments)
echo "=== fuzzy token search for 'pernicios' in col2/col10 of final DB"
awk -F'\t' '
  tolower($2) ~ /pernicios/ || (NF>=10 && tolower($10) ~ /pernicios/) {print NR "\t" $2 "\t" $10}
' "$DB" | head || true
# === fuzzy token search for 'pernicios' in col2/col10 of final DB
# 5079    Encarsia perniciosi     Encarsia perniciosi
# 10647   Phlebotomus perniciosus Phlebotomus perniciosus
# 301161  Maladera perniciosa     Maladera perniciosa
# 508280  Encarsia perniciosi     Encarsia perniciosi
# 510153  Encarsia perniciosi     Encarsia perniciosi
# 511686  Encarsia perniciosi     Aphelinidae sp. BOLD-2016
# 514683  Encarsia perniciosi     Encarsia perniciosi
# 523740  Encarsia perniciosi     Encarsia perniciosi
# 524024  Encarsia perniciosi     Aphelinidae sp. BOLD-2016
# 531320  Encarsia perniciosi     Encarsia perniciosi




##
SRC="processed/anml_amplicons_primary.txt"

# exact, case-insensitive match in col2 or col10 for all common names
for name in \
  "Diaspidiotus perniciosus" "Quadraspidiotus perniciosus" \
  "Comstockaspis perniciosa" "Aspidiotus perniciosus"
do
  echo ">> pre-pairwise ($SRC): $name"
  awk -F'\t' -v IGNORECASE=1 -v q="$name" '
    tolower($2)==tolower(q) || (NF>=10 && tolower($10)==tolower(q)) {print NR "\t" $0}
  ' "$SRC" | head || true
done

# fuzzy token to catch odd genus placements / typos
awk -F'\t' '
  tolower($2) ~ /pernicios/ || (NF>=10 && tolower($10) ~ /pernicios/) {print NR "\t" $2 "\t" $10}
' "$SRC" | head || true

