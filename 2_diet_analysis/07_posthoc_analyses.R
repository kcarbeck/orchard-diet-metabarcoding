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

# ensure factors you’ll use exist
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
# (you can still rarefy in R if a reviewer insists; otherwise use these)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
alpha_tab <- phyloseq::estimate_richness(ps, measures = c("Observed","Shannon")) %>%
  rownames_to_column("SampleID")
alpha_df <- alpha_tab %>%
  left_join(sample_data(ps) %>% as.data.frame() %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  select(SampleID, Observed, Shannon, Year, Site, Plate)
write.csv(alpha_df, file = "analysis/alpha_metrics.csv", row.names = FALSE)
ggplot(alpha_df, aes(Year, Shannon)) + geom_boxplot()

# ---- 4) beta diversity (unrarefied counts) ----
# jaccard (presence/absence) and bray-curtis
mat <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) mat <- t(mat)
dist_jac <- vegdist(mat > 0, method = "jaccard")
dist_bra <- vegdist(mat,    method = "bray")

# ordination
pcoa_bra <- cmdscale(dist_bra, k = 2, eig = TRUE)
# plot with ggplot if desired…

# PERMANOVA controlling for Plate and Site (if available)
meta <- data.frame(sample_data(ps))
if ("Site" %in% names(meta)) {
  adonis2(dist_bra ~ Year + Site + Plate, data = meta)
} else {
  adonis2(dist_bra ~ Year + Plate, data = meta)
}

# ---- 5) FOO & RRA at chosen ranks (Order=4, Family=5, Genus=6) ----
collapse_to <- function(ps, level) {
  # collapse by taxonomy rank using phyloseq
  tx <- tax_table(ps)
  rank_name <- colnames(tx)[level]
  tax_glom(ps, taxrank = rank_name, NArm = TRUE)
}

ps_Fam <- collapse_to(ps, 5)
ps_Gen <- collapse_to(ps, 6)

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

# Save Year×Site summaries (Family)
rra_group <- rra_long %>% group_by(Year, Site, Taxon) %>% summarize(RRA = mean(RRA), .groups = "drop")
write.csv(rra_group, file = "analysis/RRA_Family_YearSite.csv", row.names = FALSE)
write.csv(foo_long, file = "analysis/FOO_Family_YearSite.csv", row.names = FALSE)

# ---- 6) simple figures (customize as needed) ----
# stacked bar by sample (RRA), faceted by Year
ggplot(rra_long, aes(Sample, RRA, fill = Taxon)) + geom_col() + facet_wrap(~Year, scales="free_x")

# FOO heatmap by Year×Site (facet by Site if present)
p_foo <- ggplot(foo_long, aes(Year, Taxon, fill = FOO)) + geom_tile()
if (!all(is.na(foo_long$Site))) p_foo <- p_foo + facet_wrap(~Site, scales = "free_x")
p_foo

# ---------------- Pest-layer summaries (optional) ----------------
if (file.exists(pest_csv)) {
  pest_df <- suppressMessages(readr::read_csv(pest_csv, show_col_types = FALSE))
  if (all(c("rank","taxon") %in% names(pest_df))) {
    sanitize <- function(x) {
      x <- as.character(x)
      x <- gsub("^(k__|p__|c__|o__|f__|g__|s__)", "", x)
      trimws(tolower(x))
    }
    tx <- as.data.frame(tax_table(ps))
    tx$Species_clean <- sanitize(tx$Species)
    tx$Genus_clean   <- sanitize(tx$Genus)

    p_species <- pest_df %>% filter(tolower(rank) == "species") %>% pull(taxon) %>% sanitize()
    p_genus   <- pest_df %>% filter(tolower(rank) == "genus")   %>% pull(taxon) %>% sanitize()

    is_pest <- (tx$Species_clean %in% p_species) | (tx$Genus_clean %in% p_genus)
    names(is_pest) <- rownames(tx)

    counts <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) counts <- t(counts)
    col_tot <- pmax(1, colSums(counts))
    rra_mat <- sweep(counts, 2, col_tot, "/")
    pest_idx <- which(is_pest[rownames(rra_mat)])
    rra_pest <- if (length(pest_idx)) colSums(rra_mat[pest_idx, , drop = FALSE]) else rep(0, ncol(rra_mat))
    foo_pest <- if (length(pest_idx)) as.integer(colSums(counts[pest_idx, , drop = FALSE] > 0) > 0) else rep(0L, ncol(counts))

    smd <- sample_data(ps) %>% as.data.frame() %>% rownames_to_column("SampleID")
    pest_sample <- tibble(
      SampleID = colnames(rra_mat),
      RRA_pest = as.numeric(rra_pest),
      FOO_pest = as.integer(foo_pest)
    ) %>% left_join(smd, by = "SampleID") %>% select(SampleID, Year, Site, Plate, RRA_pest, FOO_pest)
    write.csv(pest_sample, file = "analysis/pest_sample_summary.csv", row.names = FALSE)

    pest_group <- pest_sample %>% group_by(Year, Site) %>% summarize(mean_RRA_pest = mean(RRA_pest), FOO_pest = mean(FOO_pest), .groups = "drop")
    write.csv(pest_group, file = "analysis/pest_YearSite_summary.csv", row.names = FALSE)

    # heatmap of mean_RRA_pest by Year×Site
    ggplot(pest_group, aes(Year, Site, fill = mean_RRA_pest)) + geom_tile() + scale_fill_viridis_c()

    # top pest taxa by group (Species if available else Genus)
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
}

# ---- 7) differential abundance example: ANCOM-BC2 at Family ----
otu <- as.matrix(otu_table(ps_Fam)); if (!taxa_are_rows(ps_Fam)) otu <- t(otu)
smd <- data.frame(sample_data(ps_Fam))
fit <- ancombc2(otu_table = otu, col_data = smd, fix_formula = "Year + Plate",
                p_adj_method = "BH", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05)
da_table <- fit$res %>% filter(q_val < 0.05)
write.csv(da_table, "analysis/DA_Family_ANCOMBC2.csv", row.names=FALSE)
