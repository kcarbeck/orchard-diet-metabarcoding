# evaluate a naive bayes classifier for orchard fecal samples
# katherine carbeck
# 15 oct 2025

# Species in the Woodpecker mock community (PCR +):
# Pholcus phalangioides
# Forficula auricularia
# Oecanthus niveus
# Alydus eurinus
# Vespula maculifrons
# Conocephalus brevipennis
# Eurosta solidaginis
# Sarcophaga sp.
# Xylocopa virginica

# inputs from the trained reference set 
REF_SEQS=/workdir/kcarbeck/final_outputs/ref_derep_seqs_141025.qza
REF_TAX=/workdir/kcarbeck/final_outputs/ref_derep_tax_141025.qza
CLASSIFIER=/workdir/kcarbeck/classifier/nb_anml_orchard_classifier_141025.qza
OUTDIR=/workdir/kcarbeck/classifier_out/mock_eval
mkdir -p "$OUTDIR"

# species terms to include for the Woodpecker mock community (used to filter both sequences and taxonomy
# "Sarcophaga" kept at genus level
MOCK_INCLUDE_TERMS="Pholcus_phalangioides,Forficula_auricularia,Oecanthus_niveus,Alydus_eurinus,Vespula_maculifrons,Conocephalus_brevipennis,Eurosta_solidaginis,Sarcophaga,Xylocopa_virginica"

# check for epithet (species name only)
for epithet in phalangioides auricularia niveus eurinus maculifrons brevipennis solidaginis virginica
do
  printf "%-14s  " "$epithet"; \
  awk -F'\t' 'NR>1 && tolower($2) ~ tolower("'"$epithet"'") {c++} END{print c+0}' \
    /workdir/kcarbeck/tmp_tax_export/taxonomy.tsv
done

# phalangioides   13
# auricularia     12
# niveus          9
# eurinus         7
# maculifrons     34
# brevipennis     7
# solidaginis     33
# virginica       20


##############################################################################
#*               1. subset ref seq to the mock taxa
##############################################################################
qiime taxa filter-seqs \
  --i-sequences "$REF_SEQS" \
  --i-taxonomy "$REF_TAX" \
  --p-include "$MOCK_INCLUDE_TERMS" \
  --o-filtered-sequences "$OUTDIR/ref_mock_seqs.qza"

##############################################################################
#*      2. classify those mock reference sequences with the classifier
##############################################################################
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads "$OUTDIR/ref_mock_seqs.qza" \
  --o-classification "$OUTDIR/nb_on_mock_ref.qza"

##############################################################################
#*    3. build the expected taxonomy for mock community
##############################################################################
# export full ref taxonomy
qiime tools export \
  --input-path "$REF_TAX" \
  --output-path "$OUTDIR/ref_tax_export"

# export kept mock seqs to get their feature ids
qiime tools export \
  --input-path "$OUTDIR/ref_mock_seqs.qza" \
  --output-path "$OUTDIR/ref_mock_seqs_export"

# get ids from fasta headers
grep '^>' "$OUTDIR/ref_mock_seqs_export/dna-sequences.fasta" | sed 's/^>//' > "$OUTDIR/mock_ids.txt"

# keep header + rows whose first column (feature id) is in mock_ids.txt
(head -n 1 "$OUTDIR/ref_tax_export/taxonomy.tsv" && \
  awk 'NR==FNR{ids[$1];next} (NR==1)||($1 in ids)' \
      "$OUTDIR/mock_ids.txt" "$OUTDIR/ref_tax_export/taxonomy.tsv") \
  > "$OUTDIR/expected_mock_taxonomy.tsv"

# re-import as FeatureData[Taxonomy]
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format TSVTaxonomyFormat \
  --input-path "$OUTDIR/expected_mock_taxonomy.tsv" \
  --output-path "$OUTDIR/expected_mock_tax.qza"

##############################################################################
#*    4. evaluate taxonomy: observed vs expected on the same feature ids
##############################################################################
qiime quality-control evaluate-taxonomy \
  --i-observed-taxa "$OUTDIR/nb_on_mock_ref.qza" \
  --i-expected-taxa "$OUTDIR/expected_mock_tax.qza" \
  --p-depth 7 \
  --o-visualization "$OUTDIR/eval_nb_mock_141025.qzv"

##############################################################################
#*    5. sanity checks
##############################################################################
# how many features in observed vs expected (should be equal or very close)
qiime metadata tabulate \
  --m-input-file "$OUTDIR/nb_on_mock_ref.qza" \
  --o-visualization "$OUTDIR/nb_on_mock_ref.qzv"

# count rows in expected taxonomy (minus header)
tail -n +2 "$OUTDIR/expected_mock_taxonomy.tsv" | wc -l
# 386

# output plot shows:
#ranks 1–5 (kingdom → family):
# precision, recall, and F all ≈ 1  -> classifier perfectly assigns higher-level ranks

# rank 6 (genus):
# a noticeable drop, recall ≈ 0.9 and precision still ≈ 1  -> classifier rarely mislabels genera (few false positives) but often gives genus-level “no call” or genus NA (false negatives)

# rank 7 (species):
# precision ≈ 0.95 – 1, recall ≈ 0.75  -> classifier is very conservative at species level—when it calls a species it’s almost always right, but it withholds calls for ~25 % of true species.
# the green F1 ≈ 0.85 means overall species-level accuracy is good but not perfect, which is typical for short COI- or ANML-length barcodes.

# copy out to storage:
cp /workdir/kcarbeck/classifier_out/mock_eval/eval_nb_mock_141025.qzv /lustre2/home/lc736_0001/diet/songbird_coi_database/classifier/classifier_out/
cp /workdir/kcarbeck/classifier_out/mock_eval/nb_on_mock_ref.qza /lustre2/home/lc736_0001/diet/songbird_coi_database/classifier/classifier_out/
cp /workdir/kcarbeck/classifier_out/mock_eval/expected_mock_tax.qza /lustre2/home/lc736_0001/diet/songbird_coi_database/classifier/classifier_out/