# train a naive bayes classifier for orchard fecal samples
# katherine carbeck
# 14 oct 2025

# activate qiime2 environment
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source /programs/miniconda3/bin/activate qiime2-amplicon-2024.10

##############################################################################
#*               1. train a naive bayes classifier
##############################################################################
# validate inputs
qiime tools validate /workdir/kcarbeck/final_outputs/ref_derep_seqs_141025.qza
qiime tools validate /workdir/kcarbeck/final_outputs/ref_derep_tax_141025.qza


# train a naive bayes classifier (output from super mode dereplication)
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /workdir/kcarbeck/final_outputs/ref_derep_seqs_141025.qza  \
  --i-reference-taxonomy /workdir/kcarbeck/final_outputs/ref_derep_tax_141025.qza \
  --o-classifier /workdir/kcarbeck/classifier/nb_anml_orchard_classifier_141025.qza \
  --verbose
#/programs/miniconda3/envs/qiime2-amplicon-2024.10/lib/python3.10/site-packages/q2_feature_classifier/classifier.py:106: UserWarning: The TaxonomicClassifier artifact that results from this method was trained using scikit-learn version 1.4.2. It cannot be used with other versions of scikit-learn. (While the classifier may complete successfully, the results will be unreliable.)

#Saved TaxonomicClassifier to: /workdir/kcarbeck/classifier/nb_anml_orchard_classifier_141025.qza

# validate inputs
qiime tools validate /workdir/kcarbeck/final_outputs/ref_derep_lca_seqs_141025.qza
qiime tools validate /workdir/kcarbeck/final_outputs/ref_derep_lca_tax_141025.qza

# train a naive bayes classifier (output from LCA mode dereplication)
/usr/bin/time -v \
  qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads /workdir/kcarbeck/final_outputs/ref_derep_lca_seqs_141025.qza \
    --i-reference-taxonomy /workdir/kcarbeck/final_outputs/ref_derep_lca_tax_141025.qza \
    --o-classifier /workdir/kcarbeck/classifier/nb_anml_orchard_classifier_141025_lca.qza \
    --verbose 2>&1 | tee /workdir/kcarbeck/classifier/train_nb_141025_lca.log

##############################################################################
#*               2. classify (reads already ANML-trimmed)
##############################################################################
# classify using super mode
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

# make taxa barplot -- NEED TO UPDATE METADATA FILE I THINK IT SHOULD ONLY BE FOR SAMPLES IN PLATE 1, also metadata format is wack
qiime taxa barplot \
  --i-table /workdir/kcarbeck/orchards/table_plate1.qza \
  --i-taxonomy /workdir/kcarbeck/classifier_out/nb_classified_taxonomy_141025.qza \
  --m-metadata-file /workdir/kcarbeck/orchards/fecal_sampling_full_metadata_2_CLEAN.tsv \
  --o-visualization /workdir/kcarbeck/classifier_out/barplot_plate1_before_filtering_141025.qzv

# classify using LCA mode
/usr/bin/time -v \
qiime feature-classifier classify-sklearn \
  --i-classifier /workdir/kcarbeck/classifier/nb_anml_orchard_classifier_141025_lca.qza\
  --i-reads /workdir/kcarbeck/orchards/rep-seqs_plate1.qza \
  --p-n-jobs 22 \
  --o-classification /workdir/kcarbeck/classifier_out/nb_classified_taxonomy_141025_lca.qza \
  --verbose 2>&1 | tee /workdir/kcarbeck/classifier/nb_classify_141025_lca.log


qiime metadata tabulate \
  --m-input-file /workdir/kcarbeck/classifier_out/nb_classified_taxonomy_141025_lca.qza \
  --o-visualization /workdir/kcarbeck/classifier_out/nb_classified_taxonomy_141025_lca.qzv

# make taxa barplot -- NEED TO UPDATE METADATA FILE I THINK IT SHOULD ONLY BE FOR SAMPLES IN PLATE 1, also metadata format is wack
qiime taxa barplot \
  --i-table /workdir/kcarbeck/orchards/table_plate1.qza \
  --i-taxonomy /workdir/kcarbeck/classifier_out/nb_classified_taxonomy_141025_lca.qza \
  --m-metadata-file /workdir/kcarbeck/orchards/fecal_sampling_full_metadata_2_CLEAN.tsv \
  --o-visualization /workdir/kcarbeck/classifier_out/barplot_plate1_before_filtering_141025_lca.qzv

##############################################################################


# copy out to storage:
cp /workdir/kcarbeck/classifier/nb_anml_orchard_classifier_141025.qza /lustre2/home/lc736_0001/diet/songbird_coi_database/classifier/
cp /workdir/kcarbeck/classifier/nb_anml_orchard_classifier_141025_lca.qza /lustre2/home/lc736_0001/diet/songbird_coi_database/classifier/
cp /workdir/kcarbeck/classifier_out/nb_classified_taxonomy_141025.qza /lustre2/home/lc736_0001/diet/songbird_coi_database/classifier/classifier_out/
cp /workdir/kcarbeck/classifier_out/nb_classified_taxonomy_141025_lca.qza /lustre2/home/lc736_0001/diet/songbird_coi_database/classifier/classifier_out/

# copy for tasos
cp /workdir/kcarbeck/classifier/nb_anml_orchard_classifier_141025.qza /lustre2/home/lc736_0001/orchards/classifier