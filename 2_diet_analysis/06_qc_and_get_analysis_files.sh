# QC and then branch into files for downstream analyses
# katherine carbeck
# 21 oct 2025

########################################################
#*  EDIT FILE NAMES AND PATHS 
########################################################


########################################################
#* 1) quick QC with controls
########################################################
# quick composition check includes controls
qiime taxa barplot \
  --i-table table_merged_nocontam.qza \
  --i-taxonomy taxonomy_merged.qza \
  --m-metadata-file all_plates_metadata.tsv \
  --o-visualization qc/taxa_barplot_with_controls.qzv

# export and make merged summary table
qiime feature-table summarize \
  --i-table table_merged_nocontam.qza \
  --m-sample-metadata-file all_plates_metadata.tsv \
  --o-visualization qc/table_with_controls_summary.qzv


########################################################
#* 2) make “analysis-only” ASV table with controls removed
########################################################
qiime feature-table filter-samples \
  --i-table table_merged_nocontam.qza \
  --m-metadata-file all_plates_metadata.tsv \
  --p-where "control_role = 'sample'" \
  --o-filtered-table table_analysis_only.qza

#keep table_merged_nocontam.qza for your records/QC
# use table_analysis_only.qza for everything downstream

# it should be arthropoda only already, but just in case
qiime taxa filter-table \
  --i-table table_analysis_only.qza \
  --i-taxonomy taxonomy_merged.qza \
  --p-include "Arthropoda" \
  --o-filtered-table table_analysis_only_arthropoda.qza


########################################################
#* 3)        TAXA BARPLOT
########################################################
qiime taxa barplot \
  --i-table table_analysis_only_arthropoda.qza \
  --i-taxonomy taxonomy_merged.qza \
  --m-metadata-file all_plates_metadata.tsv \
  --o-visualization taxa_barplot_final.qzv



########################################################
#* IMPORTANT NOTES REGARDING RAREFICATION 
# https://doi.org/10.1371/journal.pcbi.1003531
# 

# USE RAREFIED TABLE FOR: 
# - core metrics (observed features, Shannon, Jaccard, Bray–Curtis/UniFrac)
# - alpha diversity group tests that you want computed on an even-depth table (e.g., alpha-group-significance on rarefied table's metrics outputs)
# - any downstream analyses that require an even-depth table or eqaul library sizes
#NOTE: rarefying is generally accepted as not needed anymore for most analyses (see above), but there may be some cases where it is still needed.

# DO NOT USE RAREFIED TABLE FOR:
# - DEICODE/RPCA (uses unrarefied counts and iscompositional + sparse-friendly).
# - ANCOM-BC/ANCOM-BC2 (uses unrarefied counts with its own normalization).
# - taxa barplots --> refer relative abundance on unrarefied counts to avoid dropping samples
# - FOO (presence/absence) summaries (compute from unrarefied presence/absence; rarefying can randomly drop rare detections)
# - any step where keeping all reads/samples matters, or where the method handles unequal depth by design

########################################################

########################################################
#* 4)       (optional) alpha-rarefaction viz    
#* for depth selection only and visualization (no table)
########################################################
qiime diversity alpha-rarefaction \
  --i-table table_analysis_only_arthropoda.qza \
  --m-metadata-file all_plates_metadata.tsv \
  --p-max-depth 50000 \
  --o-visualization alpha_rarefaction_arthropoda.qzv