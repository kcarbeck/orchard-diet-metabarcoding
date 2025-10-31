

# -----------------------------
# paths
# -----------------------------
EVAL=/workdir/kcarbeck/classifier_out/mock_eval                    # or wherever your eval is
EXP_QZA="$EVAL/expected_mock_tax.qza"                    # expected taxa used in evaluate-taxonomy
OBS_QZA="$EVAL/nb_on_mock_ref.qza"                  # observed taxa output from classify-sklearn
WORK="$EVAL/confusion_out"
mkdir -p "$WORK"

# -----------------------------
# 1) export both taxonomy artifacts
#    (qiime sometimes writes taxonomy as taxonomy.tsv or metadata.tsv; handle both)
# -----------------------------
qiime tools export --input-path "$EXP_QZA" --output-path "$WORK/exp_export"
qiime tools export --input-path "$OBS_QZA" --output-path "$WORK/obs_export"

# normalize paths to the exported tsvs
EXP_TSV="$WORK/exp_export/taxonomy.tsv"           # expected
OBS_TSV="$WORK/obs_export/taxonomy.tsv"           # observed (or metadata.tsv)
[ -f "$OBS_TSV" ] || OBS_TSV="$WORK/obs_export/metadata.tsv"

# quick sanity check
head -n 3 "$EXP_TSV"
head -n 3 "$OBS_TSV"

# Feature ID      Taxon
# EU484502        k__NA;p__Arthropoda;c__Insecta;o__Diptera;f__Tephritidae;g__Eurosta;# s__Eurosta_solidaginis
# GU802453        k__NA;p__Arthropoda;c__Insecta;o__Hymenoptera;f__Apidae;g__Xylocopa;# s__Xylocopa_virginica
# Feature ID      Taxon   Confidence
# crabs_16921_Sarcophaga_subvicina        k__NA;p__Arthropoda;c__Insecta;o__Diptera;# f__Sarcophagidae;g__Sarcophaga   0.988967198534271
# crabs_34025_Sarcophaga  k__NA;p__Arthropoda;c__Insecta;o__Diptera;f__Sarcophagidae;# g__Sarcophaga;s__NA  0.9734892259752441

# -----------------------------
# 2) helper: collapse a semicolon-delimited taxonomy string to a given depth
#    keeps rank prefixes (k__/p__/.../s__), so strings from exp/obs remain comparable
# -----------------------------
collapse_to_depth='
function rstrip(s){ sub(/[[:space:]]+$/,"",s); return s }
function collapse(tx, depth,   i, n, out, a){
  n = split(tx, a, /;/)
  if (n==0) return tx
  out=""
  if (n < depth) depth=n
  for (i=1; i<=depth; i++){
    out = out ((i==1) ? "" : ";") a[i]
  }
  return rstrip(out)
}
'

# -----------------------------
# 3) build species-level (depth=7) confusion table for ids present in both sets
#    outputs: FeatureID \t expected_taxon@7 \t observed_taxon@7 \t match(0/1)
# -----------------------------
awk -F'\t' -v OFS='\t' -v DEPTH=7 '
function collapse(tx, depth,   i, n, out, a){
  n = split(tx, a, /;/)
  if (n == 0) return tx
  if (depth > n) depth = n
  out = a[1]
  for (i=2; i<=depth; i++) out = out ";" a[i]
  gsub(/[[:space:]]+$/, "", out)
  return out
}
NR==FNR { if (FNR==1) next; exp7[$1] = collapse($2, DEPTH); next }  # read expected
FNR==1 { next }                                                     # skip observed header
($1 in exp7) {
  obs7 = collapse($2, DEPTH)
  m = (obs7 == exp7[$1]) ? 1 : 0
  print $1, exp7[$1], obs7, m
}
' "$EXP_TSV" "$OBS_TSV" > "$WORK/confusion_depth7.tsv"

# quick peek
column -t -s $'\t' "$WORK/confusion_depth7.tsv" | head

# -----------------------------
# 4) aggregate to a tidy confusion summary at species level
#    expected \t observed \t count
# -----------------------------
# compact confusion matrix @ species
awk -F'\t' -v OFS='\t' '
{ if (FNR==1) next; c[$2"\t"$3]++ }
END{
  print "expected","observed","count"
  for (k in c) print k, c[k]
}
' "$WORK/confusion_depth7.tsv" | sort -k1,1 -k3,3nr > "$WORK/confusion_matrix_depth7.tsv"


# -----------------------------
# 5) per-species precision/recall at depth 7 (micro-averaged over features)
#    precision(exp)= TP / (TP + FP_as_this_species)
#    recall(exp)= TP / (TP + FN_for_this_species)
# -----------------------------

# matrix path
M="$WORK/confusion_matrix_depth7.tsv"

# run: builds per_species_scores.tsv next to the matrix file
python3 - "$M" <<'PY'
import csv, sys, os
m = sys.argv[1]
out = os.path.join(os.path.dirname(m), "per_species_scores.tsv")

exptotal = {}
obstotal = {}
tp = {}
fn = {}
fp = {}

with open(m, newline='') as f:
    r = csv.reader(f, delimiter='\t')
    next(r, None)  # skip header
    for exp, obs, cnt in r:
        cnt = int(cnt)
        exptotal[exp] = exptotal.get(exp, 0) + cnt
        obstotal[obs] = obstotal.get(obs, 0) + cnt
        if exp == obs:
            tp[exp] = tp.get(exp, 0) + cnt
        else:
            fn[exp] = fn.get(exp, 0) + cnt
            fp[obs] = fp.get(obs, 0) + cnt

rows = [("species","TP","FN","FP","Precision","Recall","F1")]
for s in exptotal:
    T  = tp.get(s, 0)
    FN = fn.get(s, 0)
    FP = fp.get(s, 0)
    P  = (T/(T+FP)) if (T+FP)>0 else 0.0
    R  = (T/(T+FN)) if (T+FN)>0 else 0.0
    F1 = (2*P*R/(P+R)) if (P+R)>0 else 0.0
    rows.append((s, T, FN, FP, f"{P:.3f}", f"{R:.3f}", f"{F1:.3f}"))

# sort by F1 desc and write
rows = [rows[0]] + sorted(rows[1:], key=lambda x: float(x[6]), reverse=True)
with open(out, "w", newline='') as w:
    csv.writer(w, delimiter='\t').writerows(rows)
print(out)
PY

# preview
column -t -s $'\t' "$(dirname "$M")/per_species_scores.tsv" | head -n 30

# /workdir/kcarbeck/classifier_out/mock_eval/confusion_out/per_species_scores.tsv
# species                                                                                                    TP   FN  FP  Precision  Recall  F1
# k__NA;p__Arthropoda;c__Arachnida;o__Araneae;f__Pholcidae;g__Pholcus;s__Pholcus_phalangioides               13   0   0   1.000      1.000   1.000
# k__NA;p__Arthropoda;c__Insecta;o__Diptera;f__Sarcophagidae;g__Sarcophaga;s__Sarcophaga_sarracenioides      1    0   0   1.000      1.000   1.000
# k__NA;p__Arthropoda;c__Insecta;o__Diptera;f__Sarcophagidae;g__Sarcophaga;s__Sarcophaga_sinuata             17   0   0   1.000      1.000   1.000
# k__NA;p__Arthropoda;c__Insecta;o__Diptera;f__Tephritidae;g__Eurosta;s__Eurosta_solidaginis                 8    0   0   1.000      1.000   1.000
# k__NA;p__Arthropoda;c__Insecta;o__Hymenoptera;f__Vespidae;g__Vespula;s__Vespula_maculifrons                34   0   0   1.000      1.000   1.000
# k__NA;p__Arthropoda;c__Insecta;o__Orthoptera;f__Gryllidae;g__Oecanthus;s__Oecanthus_niveus                 9    0   0   1.000      1.000   1.000
# k__NA;p__Arthropoda;c__Insecta;o__Orthoptera;f__Tettigoniidae;g__Conocephalus;s__Conocephalus_brevipennis  4    0   0   1.000      1.000   1.000
# k__NA;p__Arthropoda;c__Insecta;o__Hymenoptera;f__Apidae;g__Xylocopa;s__Xylocopa_virginica                  14   1   0   1.000      0.933   0.966 # one unclassified ASV
# k__NA;p__Arthropoda;c__Insecta;o__Diptera;f__Sarcophagidae;g__Sarcophaga;s__Sarcophaga_subvicina           20   6   0   1.000      0.769   0.870 # several genus-only calls
# k__NA;p__Arthropoda;c__Insecta;o__Diptera;f__Sarcophagidae;g__Sarcophaga;s__NA                             163  69  0   1.000      0.703   0.825 # inherent loss due to ambiguous species labels
# k__NA;p__Arthropoda;c__Insecta;o__Dermaptera;f__Forficulidae;g__Forficula;s__Forficula_auricularia         7    5   0   1.000      0.583   0.737 # multiple genus-level calls only
# k__NA;p__Arthropoda;c__Insecta;o__Hemiptera;f__Alydidae;g__Alydus;s__Alydus_eurinus                        4    3   0   1.000      0.571   0.727 # same
# k__NA;p__Arthropoda;c__Insecta;o__Diptera;f__Sarcophagidae;g__Sarcophaga;s__                               0    7   0   0.000      0.000   0.000 # no species-level match possible


# copy out to storage:
cp /workdir/kcarbeck/classifier_out/mock_eval/confusion_out/confusion_matrix_depth7.tsv /lustre2/home/lc736_0001/diet/songbird_coi_database/classifier/classifier_out/
cp /workdir/kcarbeck/classifier_out/mock_eval/confusion_out/per_species_scores.tsv /lustre2/home/lc736_0001/diet/songbird_coi_database/classifier/classifier_out/