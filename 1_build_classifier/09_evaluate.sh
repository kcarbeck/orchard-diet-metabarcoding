# evaluate qiime classifier
# katherine carbeck
# 11 sep 2025

#! messy code right now, needs cleaning up

# activate qiime2 environment
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source /programs/miniconda3/bin/activate qiime2-amplicon-2024.10


# set these to your actual run directory
OUTDIR="orchards/vsearch_out_20250912_141514"
CLASS="$OUTDIR/classification.qza"
SEARCH="$OUTDIR/search_results.qza"

# 1. check artifact types/sizes
qiime tools peek "$CLASS"
qiime tools peek "$SEARCH"

# 2. create .qzv you can open in https://view.qiime2.org
qiime metadata tabulate --m-input-file "$CLASS"  --o-visualization "$OUTDIR/classification.qzv" &
qiime metadata tabulate --m-input-file "$SEARCH" --o-visualization "$OUTDIR/search_results.qzv" &

# 3. rep-seqs summary (lengths, seqs)
if [[ -f orchards/rep-seqs_plate1.qza ]]; then
  qiime feature-table tabulate-seqs \
    --i-data orchards/rep-seqs_plate1.qza \
    --o-visualization "$OUTDIR/rep-seqs_plate1.qzv"
fi &


# 4. export taxonomy (feature id, taxon, confidence) for grepping counts
rm -rf "$OUTDIR/exported_taxonomy" "$OUTDIR/exported_search"
qiime tools export --input-path "$CLASS"  --output-path "$OUTDIR/exported_taxonomy" &
qiime tools export --input-path "$SEARCH" --output-path "$OUTDIR/exported_search" &

TAX_TSV="$OUTDIR/exported_taxonomy/taxonomy.tsv"
SRCH_TSV="$OUTDIR/exported_search/blast6.tsv"

# quick sanity checks
head -n 5 "$TAX_TSV"
head -n 3 "$SRCH_TSV"

# 5. compute a few summary stats - taxonomy qc counts (unassigned, confidence tiers)
TOTAL_FEATURES=$(awk -F'\t' 'NR>1{c++} END{print c+0}' "$TAX_TSV")
UNASSIGNED=$(awk -F'\t' 'NR>1 && tolower($2) ~ /unassigned/ {c++} END{print c+0}' "$TAX_TSV")
HIGH_CONF=$(awk -F'\t' 'NR>1 && $3 >= 0.90 {c++} END{print c+0}' "$TAX_TSV")
MID_CONF=$(awk -F'\t' 'NR>1 && $3 >= 0.70 && $3 < 0.90 {c++} END{print c+0}' "$TAX_TSV")
LOW_CONF=$(awk -F'\t' 'NR>1 && $3 < 0.70 {c++} END{print c+0}' "$TAX_TSV")

echo "total features: $TOTAL_FEATURES"
# total features: 19535
echo "unassigned: $UNASSIGNED"
#unassigned: 6713
echo "confidence >=0.90: $HIGH_CONF"
# confidence >=0.90: 18246
echo "confidence 0.70–0.89: $MID_CONF"
# confidence 0.70–0.89: 1289
echo "confidence <0.70: $LOW_CONF"
# confidence <0.70: 0


# assigned vs unassigned
ASSIGNED=$((TOTAL_FEATURES - UNASSIGNED))
printf "assigned: %d (%.1f%%)\nunassigned: %d (%.1f%%)\n" \
  "$ASSIGNED" "$(awk -v a="$ASSIGNED" -v t="$TOTAL_FEATURES" 'BEGIN{print 100*a/t}')" \
  "$UNASSIGNED" "$(awk -v u="$UNASSIGNED" -v t="$TOTAL_FEATURES" 'BEGIN{print 100*u/t}')"
#assigned: 12822 (65.6%)
# unassigned: 6713 (34.4%)


#do unassigned ASVs have any vsearch hits? i.e., did they fail because of no hits vs the hits didn't meet consensus/thresholds?
# list unassigned feature ids
awk -F'\t' 'NR>1 && tolower($2) ~ /unassigned/ {print $1}' "$TAX_TSV" > "$OUTDIR/unassigned_ids.txt"

# count how many of those appear in the per-hit table
# (robustly detect the query id column name)
awk -F'\t' '
  NR==FNR {u[$1]=1; next}
  NR==1 {for(i=1;i<=NF;i++) h[$i]=i; q = h["query"] ? h["query"] : (h["qseqid"] ? h["qseqid"] : (h["featureid"] ? h["featureid"] : 1)); next}
  (q && ($q in u)) {seen[$q]=1}
  END {
    for(id in u) total++
    for(id in seen) hit++
    printf "unassigned with >=1 hit rows: %d\nunassigned with 0 hit rows: %d\n", hit+0, total-(hit+0)
  }
' "$OUTDIR/unassigned_ids.txt" "$SRCH_TSV"
#unassigned with >=1 hit rows: 0
#unassigned with 0 hit rows: 6713

# features with >=1 hit overall
HIT_ANY=$(awk -F'\t' '$2!="*"{seen[$1]=1} END{for(k in seen)c++; print c+0}' "$SRCH_TSV")
printf "features with >=1 hit: %d (%.1f%%)\nfeatures with 0 hits: %d (%.1f%%)\n" \
  "$HIT_ANY" "$(awk -v h="$HIT_ANY" -v t="$TOTAL_FEATURES" 'BEGIN{print 100*h/t}')" \
  "$((TOTAL_FEATURES-HIT_ANY))" "$(awk -v h="$HIT_ANY" -v t="$TOTAL_FEATURES" 'BEGIN{print 100*(t-h)/t}')"
#features with >=1 hit: 12822 (65.6%)
#features with 0 hits: 6713 (34.4%)



# 6. taxonomy depth + species/genus label counts 
# depth by number of semicolons (rough proxy for how deep the taxonomy goes)
awk -F'\t' 'NR>1 {n=gsub(/;/,";",$2); a[n]++} END{print "semicolons\tfeatures"; for(k in a) printf "%s\t%s\n",k,a[k]}' "$TAX_TSV" | sort -n \
  > "$OUTDIR/taxonomy_depth_counts.tsv"
column -t "$OUTDIR/taxonomy_depth_counts.tsv"
# 0           6713
# semicolons  features
# 3           77
# 4           220
# 5           4069
# 6           8456



#! UPDATE: counts for s__ and genus-only (works if your labels use k__/p__/.../s__ style)
S_COUNT=$(awk -F'\t' 'NR>1 && $2 ~ /(^|; )s__/ {c++} END{print c+0}' "$TAX_TSV")
G_ONLY=$(awk -F'\t' 'NR>1 && $2 ~ /(^|; )g__/ && $2 !~ /(^|; )s__/ {c++} END{print c+0}' "$TAX_TSV")
echo "features with s__ label: $S_COUNT"
echo "genus-only labels (g__ no s__): $G_ONLY"


#! 7.scan vsearch hits for any sub-threshold matches (should be none)
# show header so we know column names
head -n 1 "$SRCH_TSV"

#! detect common identity/coverage headers and print any rows with identity<97 or cov<80
awk -F'\t' '
  NR==1{
    for(i=1;i<=NF;i++){h[$i]=i}
    pid = (h["perc_identity"] ? h["perc_identity"] : (h["pident"] ? h["pident"] : (h["percent_identity"] ? h["percent_identity"] : 0)))
    qcv = (h["query_cov"] ? h["query_cov"] : (h["qcov"] ? h["qcov"] : (h["qcovs"] ? h["qcovs"] : 0)))
    printf "# detected identity col: %s\n# detected coverage col: %s\n", pid, qcv > "/dev/stderr"
    print; next
  }
  {
    ok=1
    if(pid>0 && ($pid+0) < 97.0) ok=0
    if(qcv>0 && ($qcv+0) < 80.0) ok=0
    if(!ok) print
  }
' "$SRCH_TSV" | head -n 12


# 8. taxa barplots
# set these if you want the stacked barplots
TABLE_QZA="orchards/table_plate1.qza"
META_TSV="orchards/fecal_sampling_full_metadata_2.txt"

# quick table/metadata check (makes a .qzv so you can confirm sample ids)
qiime feature-table summarize \
  --i-table "$TABLE_QZA" \
  --m-sample-metadata-file "$META_TSV" \
  --o-visualization "$OUTDIR/table_summary.qzv" &

# make the barplot
qiime taxa barplot \
  --i-table "$TABLE_QZA" \
  --i-taxonomy "$CLASS" \
  --m-metadata-file "$META_TSV" \
  --o-visualization "$OUTDIR/taxa-barplot.qzv" &




# 9. qc summary file:
# write the headline numbers to a text file for your records
{
  echo "qc summary for: $OUTDIR"
  echo "total features: $TOTAL_FEATURES"
  echo "unassigned: $UNASSIGNED"
  echo "confidence >=0.90: $HIGH_CONF"
  echo "confidence 0.70–0.89: $MID_CONF"
  echo "confidence <0.70: $LOW_CONF"
  echo "taxonomy depth counts: $OUTDIR/taxonomy_depth_counts.tsv"
} > "$OUTDIR/qc_summary.txt"

echo "open these in https://view.qiime2.org:"
ls -1 "$OUTDIR"/*.qzv 2>/dev/null
