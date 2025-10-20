# train a naive bayes classifier for orchard fecal samples
# katherine carbeck
# 14 oct 2025

# activate qiime2 environment
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source /programs/miniconda3/bin/activate qiime2-amplicon-2024.10

# classify 
/usr/bin/time -v \
qiime feature-classifier classify-sklearn \
  --i-classifier /workdir/kcarbeck/classifier/nb_anml_orchard_classifier_141025.qza\
  --i-reads /workdir/kcarbeck/orchards/rep-seqs_plate1.qza \
  --p-n-jobs 22 \
  --o-classification /workdir/kcarbeck/classifier_out/nb_classified_taxonomy_141025.qza \
  --verbose 2>&1 | tee /workdir/kcarbeck/classifier/nb_classify_141025.log

# sanity check table & taxonomy first
qiime feature-table summarize \
  --i-table /workdir/kcarbeck/orchards/table_plate1.qza \
  --o-visualization /workdir/kcarbeck/orchards/table_plate1_summary.qzv

qiime metadata tabulate \
  --m-input-file /workdir/kcarbeck/classifier_out/nb_classified_taxonomy_141025.qza \
  --o-visualization /workdir/kcarbeck/classifier_out/nb_classified_taxonomy_141025.qzv

# make taxa barplot- update metadata file format
qiime taxa barplot \
  --i-table /workdir/kcarbeck/orchards/table_plate1.qza \
  --i-taxonomy /workdir/kcarbeck/classifier_out/nb_classified_taxonomy_141025.qza \
  --m-metadata-file /workdir/kcarbeck/orchards/fecal_sampling_full_metadata_2_CLEAN.tsv \
  --o-visualization /workdir/kcarbeck/classifier_out/barplot_plate1_before_filtering_141025.qzv