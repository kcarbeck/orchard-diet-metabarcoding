# =============================================================================
# qc and quality filtering for diet analysis
# =============================================================================
# author: katherine carbeck
# date: 28 nov 2025
#
# this script:
#   - runs qc checks on merged data
#   - filters out low-quality samples and sequences
#   - generates rarefaction curves (diagnostic - does not discard data)
#   - applies srs normalization for alpha diversity
#
# =============================================================================

# =============================================================================
# terminology - rarefaction vs rarefying (mostly for myself because I keep getting confused)
# =============================================================================
#
# rarefaction curves:
#   - a diagnostic visualization to assess sampling depth adequacy
#   - plots observed richness vs sequencing depth
#   - does NOT discard data - just helps you understand your data
#   - always generate these for qc
#
# rarefying (subsampling):
#   - a data transformation that randomly subsamples to equal depth
#   - controversial: discards valid data (mcmurdie & holmes 2014)
#   - but: still valid for alpha diversity when depth varies (weiss et al. 2017)
#   - we use srs instead (a better alternative)
#
# our approach for this study:
#   - rarefaction curves: yes (qc step in this script)
#   - alpha diversity: srs normalization (beule & karlovsky 2020)
#   - beta diversity: unrarefied (bray-curtis handles uneven depth)
#   - foo/rra: unrarefied (presence/absence less depth-sensitive)
#   - differential abundance: unrarefied (ancom-bc has internal normalization)
#
# references:
#   - mcmurdie & holmes 2014: doi:10.1371/journal.pcbi.1003531
#   - weiss et al. 2017: doi:10.1186/s40168-017-0237-y
#   - beule & karlovsky 2020: doi:10.7717/peerj.9593 (srs method)
#   - deagle et al. 2019: doi:10.1111/mec.14734 (diet metabarcoding)
# =============================================================================
#  - table_merged.qza         (combined feature table)
#  - rep-seqs_merged.qza      (combined rep sequences)
#  - taxonomy_merged.qza      (combined taxonomy)
#  - all_plates_metadata.qzv  (metadata visualization)
# =============================================================================
# configuration - edit these for your data
# =============================================================================

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source /programs/miniconda3/bin/activate qiime2-amplicon-2024.10

WORK_DIR="/lustre2/home/lc736_0001/orchards/merged"
cd "$WORK_DIR"

# input files (from previous steps)
INPUT_TABLE="$WORK_DIR/table_merged_nocontam.qza"
#TAXONOMY="$WORK_DIR/taxonomy_merged.qza"
TAXONOMY="$WORK_DIR/taxonomy_enhanced.qza"
METADATA="$WORK_DIR/all_plates_metadata.tsv"

# output directories
QC_DIR="$WORK_DIR/qc_output"
ANALYSIS_DIR="$WORK_DIR/analysis_tables"
mkdir -p "$QC_DIR" "$ANALYSIS_DIR"


# =============================================================================
# STEP 0: filter out contaminants (from decontam output)
# =============================================================================
qiime feature-table filter-features \
    --i-table $WORK_DIR/table_merged.qza \
    --m-metadata-file $WORK_DIR/decontam_output/contaminant_feature_ids_thr_010.txt \
    --p-exclude-ids \
    --o-filtered-table $WORK_DIR/table_merged_nocontam.qza
#Saved FeatureTable[Frequency] to: /lustre2/home/lc736_0001/orchards/merged/table_merged_nocontam.qza
INPUT_TABLE="$WORK_DIR/table_merged_nocontam.qza"

# =============================================================================
# STEP 1: BLAST UNCLASSIFIED ASVs and import enhanced taxonomy
# =============================================================================
# go to helper_scripts/BLAST/blast_unclassified.sh
# then come back here to import enhanced taxonomy into qiime

# import enhanced taxonomy to QIIME2
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format TSVTaxonomyFormat \
  --input-path $WORK_DIR/blast_unclassified/enhanced_taxonomy/taxonomy_enhanced.tsv \
  --output-path $WORK_DIR/taxonomy_enhanced.qza


#define taxonomy variable
TAXONOMY="$WORK_DIR/taxonomy_enhanced.qza"

# =============================================================================
# STEP 2: initial QC visualizations (with controls still in)
# =============================================================================

# taxa barplot - useful for spotting any remaining contamination patterns
qiime taxa barplot \
    --i-table "$INPUT_TABLE" \
    --i-taxonomy "$TAXONOMY" \
    --m-metadata-file "$METADATA" \
    --o-visualization "$QC_DIR/taxa_barplot_with_controls.qzv" &
#Saved Visualization to: /lustre2/home/lc736_0001/orchards/merged/qc_output/taxa_barplot_with_controls.qzv

# table summary - check read depths across all samples including controls
qiime feature-table summarize \
    --i-table "$INPUT_TABLE" \
    --m-sample-metadata-file "$METADATA" \
    --o-visualization "$QC_DIR/table_with_controls_summary.qzv" &
#Saved Visualization to: /lustre2/home/lc736_0001/orchards/merged/qc_output/table_with_controls_summary.qzv

# >>> download and view these before continuing!!!! <<<
# check: do controls look clean? any obvious contamination patterns?
#looks good

# =============================================================================
# STEP 3: remove controls (keep real samples only)
# =============================================================================

qiime feature-table filter-samples \
    --i-table "$INPUT_TABLE" \
    --m-metadata-file "$METADATA" \
    --p-where "control_role = 'sample'" \
    --o-filtered-table "$ANALYSIS_DIR/table_samples_only.qza"
#Saved FeatureTable[Frequency] to: /lustre2/home/lc736_0001/orchards/merged/analysis_tables/table_samples_only.qza

# =============================================================================
# STEP 4: filter to arthropoda only (with class-level classification)
# =============================================================================
# removes host DNA, bacteria, fungi, etc. - keeps only target diet taxa
#
# WHY "Arthropoda;c_" instead of just "Arthropoda"?
#   - sequences classified only to phylum (Arthropoda) without class-level 
#     assignment are uninformative for diet analysis
#   - they inflate read counts without providing ecological information
#   - the ";c_" ensures we keep only sequences with at least class-level ID
#   - this is standard practice in diet metabarcoding (Deagle et al. 2019)


#! do we want Annelida? and maybe Mollusca?
qiime taxa filter-table \
    --i-table "$ANALYSIS_DIR/table_samples_only.qza" \
    --i-taxonomy "$TAXONOMY" \
    --p-include "Arthropoda;c_" \
    --o-filtered-table "$ANALYSIS_DIR/table_arthropoda_unfiltered.qza"
#Saved FeatureTable[Frequency] to: /lustre2/home/lc736_0001/orchards/merged/analysis_tables/table_arthropoda_unfiltered.qza

# summary of unfiltered arthropod table
qiime feature-table summarize \
    --i-table "$ANALYSIS_DIR/table_arthropoda_unfiltered.qza" \
    --m-sample-metadata-file "$METADATA" \
    --o-visualization "$ANALYSIS_DIR/table_arthropoda_unfiltered_summary.qzv"

# >>> view this to see sample depth distribution before filtering <<<
#some samples have very low depth

# =============================================================================
# STEP 5: collapse taxonomy to species level for rarefaction curves
# =============================================================================
# WHY collapse before rarefaction?
#   - ASV-level rarefaction inflates diversity because multiple ASVs can 
#     represent the same species (intraspecific variation, sequencing artifacts)
#   - for ecological interpretation ("have I captured all prey species diversity?"),
#     species-level curves are more meaningful
#   - this is standard practice in diet metabarcoding (Deagle et al. 2019)
#
# NOTE: we collapse to species (level 7) for rarefaction curves only
#       downstream analyses use ASV-level data with taxonomic collapse as needed

qiime taxa collapse \
    --i-table "$ANALYSIS_DIR/table_arthropoda_unfiltered.qza" \
    --i-taxonomy "$TAXONOMY" \
    --p-level 7 \
    --o-collapsed-table "$ANALYSIS_DIR/table_arthropoda_species_collapsed.qza"
# Saved FeatureTable[Frequency] to: /lustre2/home/lc736_0001/orchards/merged/analysis_tables/table_arthropoda_species_collapsed.qza

# =============================================================================
# STEP 6: rarefaction curves (BEFORE setting thresholds)
# =============================================================================
# this helps you decide MIN_SAMPLE_READS and SRS_DEPTH
# look for where curves plateau - that's your target depth

# NOTE: we run rarefaction on the SPECIES-collapsed table (not ASV table)
#       this shows true species-level diversity accumulation

qiime diversity alpha-rarefaction \
    --i-table "$ANALYSIS_DIR/table_arthropoda_species_collapsed.qza" \
    --m-metadata-file "$METADATA" \
    --p-min-depth 100 \
    --p-max-depth 50000 \
    --p-steps 20 \
    --o-visualization "$QC_DIR/rarefaction_curves_species.qzv"
# Saved Visualization to: /lustre2/home/lc736_0001/orchards/merged/qc_output/rarefaction_curves_species.qzv

#! >>> STOP HERE AND REVIEW <<<
# download rarefaction_curves_species.qzv and view in qiime2 viewer
# look for:
#   - where does diversity plateau? --> this informs SRS_DEPTH
#   - how many samples have very low depth? --> this informs MIN_SAMPLE_READS
#   - if curves still rising at max depth, your data may be undersampled


# =============================================================================
# STEP 7: SET THRESHOLDS (after reviewing rarefaction curves)
# =============================================================================
# adjust these based on what you saw in the rarefaction curves

# minimum reads per sample - samples below this are dropped (set where rarefaction curves plateau)
MIN_SAMPLE_READS=10000 

# minimum total reads per ASV across all samples; removes very rare sequences that may be errors (usually pcr errors)
MIN_ASV_READS=10

# minimum number of samples an ASV must appear in
# singletons (ASV in only 1 sample) are often spurious
# set to 1 to keep all, 2+ to be more stringent (keeping all for now)
MIN_ASV_SAMPLES=1 

# SRS target depth for alpha diversity normalization -- set to where most shannon diversity index values plateau in rarefaction curve (use 15k for now)
SRS_DEPTH=15000 


# =============================================================================
# STEP 8: apply ASV filters
# =============================================================================

# filter out rare ASVs with very few total reads
# MIN_ASV_READS=10
qiime feature-table filter-features \
    --i-table "$ANALYSIS_DIR/table_arthropoda_unfiltered.qza" \
    --p-min-frequency "$MIN_ASV_READS" \
    --o-filtered-table "$ANALYSIS_DIR/table_arthropoda_freq_filtered.qza"
#Saved FeatureTable[Frequency] to: /lustre2/home/lc736_0001/orchards/merged/analysis_tables/table_arthropoda_freq_filtered.qza

# filter out singletons (ASVs in only 1 sample)
# qiime feature-table filter-features \
#    --i-table "$ANALYSIS_DIR/table_arthropoda_freq_filtered.qza" \
#    --p-min-samples "$MIN_ASV_SAMPLES" \
#    --o-filtered-table "$ANALYSIS_DIR/table_arthropoda_asv_filtered.qza"


# =============================================================================
# STEP 9: apply sample depth filter
# =============================================================================

qiime feature-table filter-samples \
    --i-table "$ANALYSIS_DIR/table_arthropoda_freq_filtered.qza" \
    --p-min-frequency "$MIN_SAMPLE_READS" \
    --o-filtered-table "$ANALYSIS_DIR/table_arthropoda.qza"
#Saved FeatureTable[Frequency] to: /lustre2/home/lc736_0001/orchards/merged/analysis_tables/table_arthropoda.qza

# summary of final filtered table
qiime feature-table summarize \
    --i-table "$ANALYSIS_DIR/table_arthropoda.qza" \
    --m-sample-metadata-file "$METADATA" \
    --o-visualization "$ANALYSIS_DIR/table_arthropoda_summary.qzv"
#Saved Visualization to: /lustre2/home/lc736_0001/orchards/merged/analysis_tables/table_arthropoda_summary.qzv
#! LOOK AT THESE OUTPUTS WHEN I GET BACK TO THE COMPUTER!

# >>> view this to check how many samples/ASVs remain <<<
# if too many samples dropped, lower MIN_SAMPLE_READS and re-run from step 8


# =============================================================================
# STEP 10: final taxa barplot
# =============================================================================

qiime taxa barplot \
    --i-table "$ANALYSIS_DIR/table_arthropoda.qza" \
    --i-taxonomy "$WORK_DIR/taxonomy_enhanced.qza" \
    --m-metadata-file "$METADATA" \
    --o-visualization "$ANALYSIS_DIR/taxa_barplot_final.qzv"
#Saved Visualization to: /lustre2/home/lc736_0001/orchards/merged/analysis_tables/taxa_barplot_final.qzv

# =============================================================================
# NOTE: SRS normalization is done in R (03_posthoc_analyses/)
# =============================================================================
# The qiime2 SRS plugin wasn't available, so we do SRS in R instead.
# SRS is only needed for alpha diversity - all other analyses use unrarefied data.
# =============================================================================





# =============================================================================
# OUTPUT SUMMARY
# =============================================================================
#
# QC outputs (in qc_output/):
#   taxa_barplot_with_controls.qzv      - check for contamination patterns
#   table_with_controls_summary.qzv     - sample depths with controls
#   rarefaction_curves_species.qzv      - diagnostic: species-level diversity curves
#
# Analysis tables (in analysis_tables/):
#   table_arthropoda.qza                    - final filtered, unrarefied (ASV-level)
#   table_arthropoda_species_collapsed.qza  - species-collapsed (for rarefaction only)
#   taxa_barplot_final.qzv                  - final composition visualization
#
# Filters applied:
#   1. removed contaminants (from decontam)
#   2. removed controls (negatives and mocks)
#   3. kept only Arthropoda sequences WITH class-level classification ("Arthropoda;c_")
#   4. removed ASVs with < MIN_ASV_READS total reads
#   5. removed samples with < MIN_SAMPLE_READS reads
#
# Next: proceed to 03_posthoc_analyses/
#   - SRS normalization is done there 
#   - Alpha diversity at species and genus levels
#   - Beta diversity, FOO/RRA all calculated there
# =============================================================================


