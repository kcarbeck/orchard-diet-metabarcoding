#katherine carbeck
# 21 oct 2025
# posthoc diet analyses

library(qiime2R)
library(phyloseq)
library(tidyverse)
library(vegan)
library(ANCOMBC) # install if needed
library(tibble)
########################################################
# THIS SCRIPT:
# - imports the ASV table and taxonomy directly from qza files
# - performs quick library-size QC
# - calculates alpha metrics
# - calculates beta diversity metrics
# - performs FOO and RRA analyses
# - performs differential abundance analysis
########################################################

# ensure output directory exists
dir.create("analysis", showWarnings = FALSE, recursive = TRUE)

# optional: path to a CSV of orchard pests. Expected columns:
#   rank ("Species" or "Genus"), taxon (e.g., "Cydia pomonella" or "Rhagoletis")
# leave as-is if you don't want pest-layer summaries
pest_csv <- "data/pest_taxa.csv"

# ---- 1) import straight from qza ----
ps <- qza_to_phyloseq(
  features = "table_analysis_only_arthropoda.qza",
  taxonomy = "taxonomy_merged.qza",
  metadata = "all_plates_metadata.tsv"
)

# create df of sample read counts by sample (ensure factors you’ll use exist)
#! might need coerce to data.frame?
sample_data(ps)$Year  <- as.factor(sample_data(ps)$Year)
sample_data(ps)$Plate <- as.factor(sample_data(ps)$Plate)
if (!is.null(sample_data(ps)$Site)) {
  sample_data(ps)$Site  <- as.factor(sample_data(ps)$Site)
}

# ---- 2) quick library-size QC (samples only) ----
lib <- data.frame(SampleID = sample_names(ps),
                  reads = sample_sums(ps),
                  Year = sample_data(ps)$Year,
                  Plate = sample_data(ps)$Plate)

# ---- 3) alpha metrics (no rarefaction) ----
# how varied are each samples diet?
#! NOTES TO SELF: 
#! MAY BE WISE TO COLLAPSE TO SPECIES/GENUS LEVEL FIRST FOR ALPHA DIVRESITY (OBSERVED COUNT): ASV richness might slightly overestimate unique prey taxa. collapse ASVs to species or genus level and then count unique taxa per sample.
# (you can still rarefy in R if a reviewer insists; otherwise use these)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
alpha_tab <- phyloseq::estimate_richness(ps, measures = c("Observed","Shannon")) %>%
  rownames_to_column("SampleID")
alpha_df <- alpha_tab %>%
  left_join(sample_data(ps) %>% as.data.frame() %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  select(SampleID, Observed, Shannon, Year, Site, Plate)
write.csv(alpha_df, file = "analysis/alpha_metrics.csv", row.names = FALSE)
ggplot(alpha_df, aes(Site, Shannon)) + geom_boxplot()
ggplot(alpha_df, aes(Year, Shannon)) + geom_boxplot()

# add: non-parametric test or linear model to test if Year or Site affects alpha diversity 
# wilcoxon test for shannon by year or site?

# This approach is somewhat debated: richness is sensitive to depth, so non-rarefied observed ASVs will correlate with sample read count. Shannon is less sensitive but can still be mildly influenced by depth. Many ecologists still prefer to rarefy for alpha diversity to be safe, or use an estimated richness (like Chao1) that tries to account for unseen diversity. The pipeline’s stance is to avoid unnecessary data loss, which is a more modern viewpoint. To mitigate depth bias, they could incorporate an offset or model if comparing richness across groups (e.g., use read count as a covariate in an ANOVA). The pipeline doesn’t go that far, but it does allow the user to check rarefaction curves to ensure most samples are near-saturated (if all samples had, say, >5,000 reads and the rarefaction curves plateaued by 2,000, then even the lowest-depth sample might capture most of its diversity, making non-rarefied comparisons acceptable).


# ---- 4) beta diversity (unrarefied counts) ----
# how do diets differ among groups?
# jaccard (presence/absence) - emphasize differences in prey composition
mat <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) mat <- t(mat)
dist_jac <- vegdist(mat > 0, method = "jaccard")

#bray-curtis (abundance-based) - incorporates how dominant each prey type was
dist_bra <- vegdist(mat, method = "bray")

# ordination
pcoa_bra <- cmdscale(dist_bra, k = 2, eig = TRUE)
# plot with ggplot if desired…

# PERMANOVA on bray-curtis distance - use for controlling for plate, site, and/or year as factors. allows us to test for effects of year,for example, by partioning out the effects of plate and site
# consider an interaction (year x site) if looking for different yearly trends in different orchards?
# 
meta <- data.frame(sample_data(ps))
if ("Site" %in% names(meta)) {
  adonis2(dist_bra ~ Year + Site + Plate, data = meta)
} else {
  adonis2(dist_bra ~ Year + Plate, data = meta)
}

# ---- 5) Frequency of occurence (FOO) & Relative read abundance (RRA) --------------
# taxonomic summaries at chosen ranks (Order=4, Family=5, Genus=6, Species=7) ---> what are the main prey taxa and their frequencies?
# FOO = fraction of samples (many times within a group, like Year or Site) in which a given taxon occurs 
# RRA = proportion of sequence reads in a sample (or group of samples) that belong to a given taxon 
# NOTE: both these metrics have pros and cons, and it’s considered best practice to examine both in metabarcoding diets (Deagle et al. 2019: recommend caution in over-interpreting read counts and suggest using occurrence metrics as a complement) https://doi.org/10.1111/mec.14734

#! NOTE TO SELF: add overall FOO across all samples (not just by group) to see which prey are most universally found vs whos rare

collapse_to <- function(ps, level) {
  # collapse by taxonomy rank using phyloseq
  tx <- tax_table(ps)
  rank_name <- colnames(tx)[level]
  tax_glom(ps, taxrank = rank_name, NArm = TRUE)
}

# ps_Ord <- collapse_to(ps, 4)
ps_Fam <- collapse_to(ps, 5)
ps_Gen <- collapse_to(ps, 6)
# ps_Spe <- collapse_to(ps, 7)

# RRA per sample
ps_Fam_rel <- transform_sample_counts(ps_Fam, function(x) x / sum(x))
tmp_rra <- psmelt(ps_Fam_rel)
if (!"Site" %in% names(tmp_rra)) tmp_rra$Site <- NA
rra_long <- tmp_rra %>% transmute(Sample=Sample, Taxon=Family, RRA=Abundance, Year, Site)

# FOO per group (mean presence)
tmp_foo <- psmelt(ps_Fam)
if (!"Site" %in% names(tmp_foo)) tmp_foo$Site <- NA
foo_long <- tmp_foo %>%
  transmute(Sample=Sample, Taxon=Family, Present = as.integer(Abundance>0), Year, Site) %>%
  group_by(Year, Site, Taxon) %>%
  summarize(FOO = mean(Present), .groups="drop")

# Save Year×Site summaries (FAMILY LEVEL)
rra_group <- rra_long %>% group_by(Year, Site, Taxon) %>% summarize(RRA = mean(RRA), .groups = "drop")
write.csv(rra_group, file = "analysis/RRA_Family_YearSite.csv", row.names = FALSE)
write.csv(foo_long, file = "analysis/FOO_Family_YearSite.csv", row.names = FALSE)
# Example: Coleoptera occurred in 50% of fecal samples and constituted on average 30% of the reads in year 2025, whereas in 2024 it occurred in 80% of samples…” etc

# ---- 6) simple figures (customize as needed) ----
# stacked bar by sample (RRA), faceted by Year
ggplot(rra_long, aes(Sample, RRA, fill = Taxon)) + geom_col() + facet_wrap(~Year, scales="free_x")

# FOO heatmap by Year×Site (facet by Site if present)
p_foo <- ggplot(foo_long, aes(Year, Taxon, fill = FOO)) + geom_tile()
if (!all(is.na(foo_long$Site))) p_foo <- p_foo + facet_wrap(~Site, scales = "free_x")
p_foo

# ---------------- Pest-layer summaries ----------------
# are birds consuming known pests? how prevalent are pests in their diet?
pest_df <- readr::read_csv(pest_csv, show_col_types = FALSE)
if (all(c("rank","taxon") %in% names(pest_df))) {
  sanitize <- function(x) {
    x <- as.character(x)
    x <- gsub("^(k__|p__|c__|o__|f__|g__|s__)", "", x)
    trimws(tolower(x))
  }
  tx <- as.data.frame(tax_table(ps))
  tx$Species_clean <- sanitize(tx$Species)
  tx$Genus_clean   <- sanitize(tx$Genus)  p_species <- pest_df %>% filter(tolower(rank) == "species") %>% pull(taxon) %>% sanitize()
  p_genus   <- pest_df %>% filter(tolower(rank) == "genus")   %>% pull(taxon) %>% sanitize()  is_pest <- (tx$Species_clean %in% p_species) | (tx$Genus_clean %in% p_genus)
  names(is_pest) <- rownames(tx)  counts <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) counts <- t(counts)
  col_tot <- pmax(1, colSums(counts))
  rra_mat <- sweep(counts, 2, col_tot, "/")
  pest_idx <- which(is_pest[rownames(rra_mat)])
  rra_pest <- if (length(pest_idx)) colSums(rra_mat[pest_idx, , drop = FALSE]) else rep(0, ncol(rra_mat))
  foo_pest <- if (length(pest_idx)) as.integer(colSums(counts[pest_idx, , drop = FALSE] > 0) > 0) else rep(0L, ncol(counts))  smd <- sample_data(ps) %>% as.data.frame() %>% rownames_to_column("SampleID")
  pest_sample <- tibble(
    SampleID = colnames(rra_mat),
    RRA_pest = as.numeric(rra_pest),
    FOO_pest = as.integer(foo_pest)
  ) %>% left_join(smd, by = "SampleID") %>% select(SampleID, Year, Site, Plate, RRA_pest, FOO_pest)
  write.csv(pest_sample, file = "analysis/pest_sample_summary.csv", row.names = FALSE)  pest_group <- pest_sample %>% group_by(Year, Site) %>% summarize(mean_RRA_pest = mean(RRA_pest), FOO_pest = mean(FOO_pest), .groups = "drop")
  write.csv(pest_group, file = "analysis/pest_YearSite_summary.csv", row.names = FALSE)  # heatmap of mean_RRA_pest by Year×Site
  ggplot(pest_group, aes(Year, Site, fill = mean_RRA_pest)) + geom_tile() + scale_fill_viridis_c()  # top pest taxa by group (Species if available else Genus)
  tx$feature_id <- rownames(tx)
  tx$taxon_best <- ifelse(!is.na(tx$Species_clean) & nzchar(tx$Species_clean), tx$Species_clean, tx$Genus_clean)
  if (length(pest_idx)) {
    rra_long_all <- as_tibble(rra_mat, rownames = "feature_id") %>%
      filter(feature_id %in% rownames(tx)[is_pest]) %>%
      pivot_longer(-feature_id, names_to = "SampleID", values_to = "RRA") %>%
      left_join(tx[, c("feature_id","taxon_best")], by = "feature_id") %>%
      left_join(smd, by = "SampleID")
    top_pest <- rra_long_all %>% group_by(Year, Site, taxon_best) %>% summarize(mean_RRA = mean(RRA), .groups = "drop") %>% arrange(Year, Site, desc(mean_RRA))
    write.csv(top_pest, file = "analysis/pest_top_taxa_YearSite.csv", row.names = FALSE)
  }
}

# ---- 7) differential abundance example: ANCOM-BC2 at FAMILY LEVEL ----
# test for differential abundance of prey families between groups (year, controlling for plate and/or site). E.g., use for statistical testing: "family X was significantly more abundant in diets in 2025 than 2024”
# ANCOM-BC2 has improved on old methods because it deals with overdispersion and zero-inflation better. by default it includes normalization and accounts for varying sequencing depths.
# prv_cut = 0.10 (only families present in at least 10% of samples are tested, to avoid very rare taxa) 
# lib_cut = 1000 (excluding samples with <1000 reads from the DA analysis to avoid extreme low-depth noise
otu <- as.matrix(otu_table(ps_Fam)); if (!taxa_are_rows(ps_Fam)) otu <- t(otu)
smd <- data.frame(sample_data(ps_Fam))
fit <- ancombc2(otu_table = otu, col_data = smd, fix_formula = "Year + Plate",
                p_adj_method = "BH", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05)
da_table <- fit$res %>% filter(q_val < 0.05)
write.csv(da_table, "analysis/DA_Family_ANCOMBC2.csv", row.names=FALSE)

# INTERPRET WITH CAUTION NOTE: a statistically significant increase in one prey family could mean real diet shift or could be influenced by detection biases, but given that this method tries to account for those, it’s a strong approach
