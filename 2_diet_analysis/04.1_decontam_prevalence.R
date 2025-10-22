# decontam prevalence method
# katherine carbeck
# 20 oct 2025

## TO RUN ON BIOHPC:
# type :"R" in the terminal to start R
# then can copy and paste the code below to run the script

# https://github.com/jbisanz/qiime2R/

#only install once
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # this takes a long time to compile
BiocManager::install("decontam")
library(qiime2R)     # read_qza, qza_to_phyloseq, read_q2metadata
library(phyloseq)
library(decontam)
library(ggplot2)

# =============================================================================
# Decontam (prevalence method)
# -----------------------------------------------------------------------------
# this script is designed to use the files resulting from the classifier (using 2025 plate 1 as an example)
#   - table_plate1.qza                            (feature table)
#   - nb_classified_taxonomy_141025.qza           (taxonomy annotations)
#   - fecal_sampling_metadata_plate1.txt          (sampling metadata)
#
# What this script does:
#   1) loads the feature table (counts), taxonomy, and sample metadata
#   2) builds a phyloseq object for one plate
#   3) identifies contaminants using the prevalence method in decontam
#           negatives = eblank, pblank, blank, empty
#           positives = pos (mock); kept for qc but excluded from the decontam fit
#   4) writes a list of contaminant feature IDs and a summary CSV

## === edit these file paths =================================================
# qiime2 artifacts (.qza)
feature_table_qza <- "/workdir/kcarbeck/2025/table_plate1.qza"
taxonomy_qza <- "/workdir/kcarbeck/2025/nb_classified_taxonomy_141025.qza"

# metadata (qiime2-style tsv) with columns: SampleID (or #SampleID) and a Species* column
metadata_path <- "/workdir/kcarbeck/2025/fecal_sampling_metadata_plate1.txt"

# output folder 
output_dir <- "/workdir/kcarbeck/2025/decontam_plate1"

# control labels as they appear in the metadata "species" column (case/spacing ignored)
negative_labels_in_species <- c("EBLANK","PBLANK","BLANK","EMPTY")
positive_labels_in_species <- c("POS")

# if your "species" column header varies, adjust this regex (we match the first column starting with "species")
species_col_regex <- "^\\s*Species"

# prevalence thresholds to sweep (documented and saved before you pick one)
thr_grid <- c(0.05, 0.10, 0.20, 0.30, 0.50)

# neg-only rule: flag features present in >= this many negatives and 0 real samples
# rationale: with few negatives, a feature in just 1 negative but 0 true samples is strong evidence of contaminant
min_neg_presence <- 1

# ==============================================================================
# helpers and metadata cleanup
# ==============================================================================
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# helper: normalize labels (uppercase; drop spaces/underscores/dashes/etc)
normlab <- function(x) toupper(gsub("[^A-Za-z0-9]+", "", trimws(as.character(x))))

# collapse newline characters that appear inside quoted cells in a txt file
fix_newlines_inside_quotes <- function(in_path) {
  txt <- readChar(in_path, file.info(in_path)$size, useBytes = TRUE)
  txt <- gsub("\r\n?", "\n", txt, perl = TRUE)
  ch <- strsplit(txt, "", fixed = TRUE)[[1]]
  out <- ch; inq <- FALSE
  for (i in seq_along(ch)) {
    if (ch[i] == "\"") inq <- !inq
    if (ch[i] == "\n" && inq) out[i] <- " "
  }
  tmp <- tempfile(fileext = ".tsv")
  writeChar(paste(out, collapse = ""), tmp, eos = NULL, useBytes = TRUE)
  tmp
}

# normalize sample ids (trim + collapse whitespace to underscores)
.normid <- function(x) gsub("\\s+", "_", trimws(as.character(x)))

prep_metadata_for_qza <- function(metadata_path) {
  # 1) fix embedded newlines and read with qiime2r
  meta_fixed1 <- fix_newlines_inside_quotes(metadata_path)
  md <- qiime2R::read_q2metadata(meta_fixed1)  # takes a path (not a data.frame)

  # 2) ensure we have SampleID (some sheets still use #SampleID)
  if (!("SampleID" %in% names(md)) && "#SampleID" %in% names(md)) {
    names(md)[names(md) == "#SampleID"] <- "SampleID"
  }
  if (!("SampleID" %in% names(md))) stop("metadata must include 'SampleID' (or '#SampleID')")

  # 3) clean ids: normalize, drop blanks, drop duplicates (keep first)
  md$SampleID <- .normid(md$SampleID)
  md <- md[!is.na(md$SampleID) & md$SampleID != "", , drop = FALSE]
  if (anyDuplicated(md$SampleID)) {
    dups <- unique(md$SampleID[duplicated(md$SampleID)])
    message("warning: duplicate SampleID(s) removed (kept first): ", paste(dups, collapse = ", "))
    md <- md[!duplicated(md$SampleID), , drop = FALSE]
  }

  # 4) write a cleaned copy for qza_to_phyloseq (expects a file path)
  meta_clean_path <- tempfile(fileext = ".tsv")
  write.table(md, meta_clean_path, sep = "\t", quote = FALSE, row.names = FALSE)
  meta_clean_path
}

# human-readable file tag for a threshold (e.g., "0.10" -> "010")
thr_tag_fn <- function(x) gsub("\\.", "", sprintf("%.2f", x))

# =============================================================================
# load data → build phyloseq → tag controls
# =============================================================================
meta_clean_path <- prep_metadata_for_qza(metadata_path)

ps <- qiime2R::qza_to_phyloseq(
  features = feature_table_qza,
  taxonomy = taxonomy_qza,
  metadata = meta_clean_path
)

# map controls from a "species" column in metadata
species_col <- grep(species_col_regex, colnames(phyloseq::sample_data(ps)), ignore.case = TRUE, value = TRUE)
if (!length(species_col)) stop("could not find a 'species' column in metadata (adjust species_col_regex)")

labs      <- normlab(phyloseq::sample_data(ps)[[species_col[1]]])
neg_labels <- normlab(negative_labels_in_species)
pos_labels <- normlab(positive_labels_in_species)

is_neg <- labs %in% neg_labels
is_pos <- labs %in% pos_labels

control_role <- ifelse(is_pos, "pos", ifelse(is_neg, "neg", "sample"))

sd <- phyloseq::sample_data(ps)
sd$control_role <- factor(control_role, levels = c("sample","neg","pos"))
sd$is_neg <- sd$control_role == "neg"
phyloseq::sample_data(ps) <- sd

# save a control-role map for later provenance
ctrl_map <- data.frame(SampleID = rownames(sd), control_role = as.character(sd$control_role))
write.csv(ctrl_map, file = file.path(output_dir, "sample_control_roles.csv"), row.names = FALSE)

message("control counts:")
print(table(sd$control_role, useNA = "ifany"))
if (sum(sd$is_neg) == 0) warning("no negatives detected; decontam prevalence will not be informative")
#sample    neg    pos
#    80     14      2

# exclude positives (mocks) from the model fit, per decontam guidance
ps_fit <- phyloseq::prune_samples(sd$control_role != "pos", ps)

# drop zero-read samples from the fit (decontam ignores these, but we record them)
zs <- sample_sums(ps_fit)
dropped_ids <- names(zs)[zs == 0]
if (length(dropped_ids)) {
  write.table(data.frame(SampleID = dropped_ids),
              file = file.path(output_dir, "dropped_zero_read_samples.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
}

# -----------------------------------------------------------------------------
# construct binary presence/absence per feature in negatives vs samples
# rationale: used for (a) neg-only rule and (b) interpretable prevalence plots
# -----------------------------------------------------------------------------
counts <- as(otu_table(ps_fit), "matrix")
if (!taxa_are_rows(ps_fit)) counts <- t(counts)  # ensure taxa are rows
counts_bin <- counts > 0

neg_flag <- sample_data(ps_fit)$is_neg
n_neg <- sum(neg_flag)
n_smp <- sum(!neg_flag)

neg_prev <- rowMeans(counts_bin[,  neg_flag, drop = FALSE])
smp_prev <- rowMeans(counts_bin[, !neg_flag, drop = FALSE])

# neg-only rule: present in at least min_neg_presence negatives and 0 real samples
neg_prev_counts <- rowSums(counts_bin[,  neg_flag, drop = FALSE])
smp_prev_counts <- rowSums(counts_bin[, !neg_flag, drop = FALSE])

neg_only_ids <- names(which((neg_prev_counts >= min_neg_presence) & (smp_prev_counts == 0)))


# =============================================================================
# sweep prevalence thresholds first (to inform choice)
# =============================================================================

impact_for_threshold <- function(t) {
  ct <- decontam::isContaminant(ps_fit, method = "prevalence", neg = "is_neg", threshold = t)
  ids_ct    <- rownames(ct)[ct$contaminant]
  ids_union <- union(ids_ct, neg_only_ids)

  # per-threshold prevalence bookkeeping
  in_samples <- sum((rownames(counts) %in% ids_union) & (smp_prev > 0))
  neg_only   <- sum((rownames(counts) %in% ids_union) & (smp_prev == 0))
  in_both    <- sum((rownames(counts) %in% ids_union) & (smp_prev > 0) & (neg_prev > 0))

  ps_tmp <- prune_taxa(!(taxa_names(ps_fit) %in% ids_union), ps_fit)
  rb <- sample_sums(ps_fit); ra <- sample_sums(ps_tmp)

  data.frame(
    threshold          = t,
    n_features_removed = length(ids_union),
    pct_reads_removed  = 100 * (1 - sum(ra) / sum(rb)),
    n_flagged_in_samples = in_samples,
    n_flagged_neg_only   = neg_only,
    n_flagged_in_both    = in_both
  )
}

sweep_res <- do.call(rbind, lapply(thr_grid, impact_for_threshold))
write.csv(sweep_res, file.path(output_dir, "threshold_sweep.csv"), row.names = FALSE)
print(sweep_res)
# threshold n_features_removed pct_reads_removed n_flagged_in_samples
#     0.05                 27      0.0006921138                    0
#     0.10                 27      0.0006921138                    0
#     0.20                 51      0.9577654205                   24
#     0.30                 56      1.7083861892                   29
#     0.50                 63      1.9534751753                   36
# n_flagged_neg_only n_flagged_in_both
#                 27                 0
#                 27                 0
#                 27                24
#                 27                29
#                 27                36

# quick plot of sweep (reads removed vs threshold)
p_sweep <- ggplot(sweep_res, aes(threshold, pct_reads_removed)) +
  geom_line() + geom_point() +
  labs(x = "prevalence threshold", y = "% reads removed (union rule)",
       title = "decontam threshold sweep") +
  theme_minimal()
ggsave(file.path(output_dir, "threshold_sweep.png"), p_sweep, width = 5, height = 3.5, dpi = 300)

# =============================================================================
# choose working threshold (edit if you want stricter 0.50); 0.1 + neg-only looked like a good choice here
# because there's minimal risk of tossing true low-abundance taxa
# =============================================================================

use_threshold <- 0.10
thr_tag <- thr_tag_fn(use_threshold)

# =============================================================================
# compute final outputs for chosen threshold
# =============================================================================

ct_use  <- decontam::isContaminant(ps_fit, method = "prevalence", neg = "is_neg", threshold = use_threshold)
ids_ct  <- rownames(ct_use)[ct_use$contaminant]
contam_ids <- union(ids_ct, neg_only_ids)  # union of decontam + neg-only rule

# select p-value column robustly across decontam versions
pcol  <- intersect(c("p.prev","p.freq","p"), colnames(ct_use))
pvals <- if (length(pcol)) ct_use[[pcol[1]]] else rep(NA_real_, nrow(ct_use))

# build per-feature summary table (handy for methods + supplement)
summary_df <- data.frame(
  feature_id      = rownames(counts),
  control_prev    = neg_prev[ rownames(counts) ],
  sample_prev     = smp_prev[ rownames(counts) ],
  neg_reads_sum   = rowSums(counts[,  neg_flag, drop = FALSE])[ rownames(counts) ],
  sample_reads_sum= rowSums(counts[, !neg_flag, drop = FALSE])[ rownames(counts) ],
  p_value         = pvals[ rownames(counts) ],
  q_value         = if (all(is.na(pvals))) NA_real_ else p.adjust(pvals, "BH"),
  contaminant_thr = ct_use$contaminant[ rownames(counts) ],
  neg_only_rule   = rownames(counts) %in% neg_only_ids,
  flagged_union   = rownames(counts) %in% contam_ids,
  stringsAsFactors = FALSE
)
write.csv(summary_df,
          file = file.path(output_dir, paste0("contaminants_summary_thr_", thr_tag, ".csv")),
          row.names = FALSE)
head(summary_df)
#                                                        feature_id control_prev
# 904c18ca89ae53a9140111b7c5a79343 904c18ca89ae53a9140111b7c5a79343   0.07142857
# 40a8b6047287228a7c0d6c8bf6637b56 40a8b6047287228a7c0d6c8bf6637b56   0.00000000
# bf7c357da2a4ba4cd3b85dd72ec000f3 bf7c357da2a4ba4cd3b85dd72ec000f3   0.00000000
# 4e7d8c9075c8cdf6b34a3cf993037455 4e7d8c9075c8cdf6b34a3cf993037455   0.00000000
# 5bad6e45017119a18b6dcdca25de484e 5bad6e45017119a18b6dcdca25de484e   0.00000000
# 59ebba2d705328de40f6d2c09d7cb367 59ebba2d705328de40f6d2c09d7cb367   0.00000000
#                                  sample_prev neg_reads_sum sample_reads_sum
# 904c18ca89ae53a9140111b7c5a79343      0.0000             2                0
# 40a8b6047287228a7c0d6c8bf6637b56      0.0125             0               76
# bf7c357da2a4ba4cd3b85dd72ec000f3      0.0125             0               11
# 4e7d8c9075c8cdf6b34a3cf993037455      0.0125             0              182
# 5bad6e45017119a18b6dcdca25de484e      0.0125             0               19
# 59ebba2d705328de40f6d2c09d7cb367      0.0125             0              664
#                                  p_value q_value contaminant_thr neg_only_rule
# 904c18ca89ae53a9140111b7c5a79343      NA      NA              NA          TRUE
# 40a8b6047287228a7c0d6c8bf6637b56      NA      NA              NA         FALSE
# bf7c357da2a4ba4cd3b85dd72ec000f3      NA      NA              NA         FALSE
# 4e7d8c9075c8cdf6b34a3cf993037455      NA      NA              NA         FALSE
# 5bad6e45017119a18b6dcdca25de484e      NA      NA              NA         FALSE
# 59ebba2d705328de40f6d2c09d7cb367      NA      NA              NA         FALSE
#                                  flagged_union
# 904c18ca89ae53a9140111b7c5a79343          TRUE
# 40a8b6047287228a7c0d6c8bf6637b56         FALSE
# bf7c357da2a4ba4cd3b85dd72ec000f3         FALSE
# 4e7d8c9075c8cdf6b34a3cf993037455         FALSE
# 5bad6e45017119a18b6dcdca25de484e         FALSE
# 59ebba2d705328de40f6d2c09d7cb367         FALSE

# write contaminant id list (one per line)
write.table(data.frame(FeatureID = contam_ids),
            file = file.path(output_dir, paste0("contaminant_feature_ids_thr_", thr_tag, ".txt")),
            quote = FALSE, sep = "\t", row.names = FALSE)

# per-sample impact for the chosen threshold
ps_clean <- prune_taxa(!(taxa_names(ps_fit) %in% contam_ids), ps_fit)
reads_before <- sample_sums(ps_fit)
reads_after  <- sample_sums(ps_clean)
# turn sample_data into a plain data.frame with rownames = sample ids
sd_df <- as.data.frame(sample_data(ps_fit))        # keeps rownames = sample IDs
role_vec <- as.character(sd_df$control_role)       # coerce to plain character labels
names(role_vec) <- rownames(sd_df)                 # name by sample id
control_role_vec <- unname(role_vec[names(reads_before)])

impact_df <- data.frame(
  sample_id    = names(reads_before),
  control_role = control_role_vec,
  reads_before = as.integer(reads_before),
  reads_after  = as.integer(reads_after),
  frac_removed = pmax(0, 1 - (reads_after / pmax(1, reads_before))),
  stringsAsFactors = FALSE
)

write.csv(impact_df,
          file = file.path(output_dir, paste0("per_sample_impact_thr_", thr_tag, ".csv")),
          row.names = FALSE)
head(impact_df)
# sample_id control_role reads_before reads_after frac_removed
#       236       sample      1024139     1024139            0
#       237       sample      1055689     1055689            0
#       238       sample       793511      793511            0
#       239       sample      1035082     1035082            0
#       240       sample       572832      572832            0
#       241       sample       899336      899336            0

# save session info and lightweight rds for reproducibility
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
saveRDS(list(ps_fit = ps_fit, ps_clean = ps_clean, contam_ids = contam_ids, sweep = sweep_res),
        file = file.path(output_dir, paste0("decontam_objects_thr_", thr_tag, ".rds")))

# =============================================================================
# diagnostic figures - can look at these in iterm2 using imgcat image.png
# =============================================================================

lib_df <- data.frame( sample_id = names(reads_before), reads = as.numeric(reads_before)) 
lib_df$control_role <- as.character(sample_data(ps_fit)[lib_df$sample_id, "control_role"][,1]) 
p_lib <- ggplot(lib_df, aes(control_role, reads)) + 
    geom_violin(fill = "grey90", color = "grey40") + 
    geom_jitter(width = 0.2, alpha = 0.4, size = 0.9) + 
    scale_y_continuous(labels = scales::comma) + 
    labs(x = "control role", y = "reads (pre-filter)", 
    title = "library sizes by control role") + theme_minimal() 
ggsave(file.path(output_dir, "library_sizes_by_control_role.png"), p_lib, width = 5, height = 3.5, dpi = 300)


# 1) library sizes by control role (confirm negatives are low-depth)
ps_for_lib <- ps_fit  

lib_df <- data.frame(
  sample_id    = sample_names(ps_for_lib),
  control_role = as.character(sample_data(ps_for_lib)$control_role),
  LibrarySize  = as.numeric(sample_sums(ps_for_lib)),
  stringsAsFactors = FALSE
)

# make sure it's an ordered factor for nice x-axis grouping
lib_df$control_role <- factor(lib_df$control_role,
                              levels = c("sample","neg","pos"))
lib_df$control_role <- droplevels(lib_df$control_role)

# color palette
pal_master <- c(sample="#0072B2", neg="#D55E00", pos="#009E73")
pal <- pal_master[levels(lib_df$control_role)]

p_lib_col <- ggplot(lib_df, aes(x = control_role, y = LibrarySize, fill = control_role)) +
  geom_violin(color = "grey20", linewidth = 0.4, alpha = 0.9, trim = FALSE) +
  geom_jitter(aes(color = control_role), width = 0.15, size = 1, alpha = 0.6) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_color_manual(values = pal, guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "control role", y = "reads (pre-filter)",
       title = "library sizes by control role") +
  theme_minimal()
ggsave(file.path(output_dir, "library_sizes_by_control_role.png"),
       p_lib_col, width = 5, height = 3.5, dpi = 300)

# 1a)another way to visualize the library sizes by control role
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
p_freq_lib<-ggplot(data=df, aes(x=Index, y=LibrarySize, color=control_role)) + geom_point()
ggsave(file.path(output_dir, "frequency_library_size_plot.png"),
       p_freq_lib, width = 5, height = 5, dpi = 300)

# 2) prevalence scatter (jittered) highlighting flagged union
plot_df <- subset(summary_df, !is.na(sample_prev) & !is.na(control_prev))
plot_df$flag <- ifelse(plot_df$flagged_union, "flagged", "kept")
p_prev_jit <- ggplot(plot_df, aes(sample_prev, control_prev, color = flag)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(size = 0.6, alpha = 0.5,
             position = position_jitter(width = 1/n_smp/2, height = 1/n_neg/2)) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "sample prevalence", y = "negative prevalence",
       title = sprintf("prevalence (jittered), thr=%s", thr_tag)) +
  theme_minimal()
ggsave(file.path(output_dir, sprintf("prevalence_scatter_jitter_thr_%s.png", thr_tag)),
       p_prev_jit, width = 5, height = 5, dpi = 300)

# 3) prevalence heat/count plot (helps when many points pile at discrete levels)
p_prev_cnt <- ggplot(plot_df, aes(sample_prev, control_prev)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_bin2d(bins = 50) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "sample prevalence", y = "negative prevalence",
       title = sprintf("prevalence (2d bins), thr=%s", thr_tag)) +
  theme_minimal()
ggsave(file.path(output_dir, sprintf("prevalence_scatter_bin2d_thr_%s.png", thr_tag)),
       p_prev_cnt, width = 5, height = 5, dpi = 300)

# 4) per-sample fraction removed (should be near-zero here; flag outliers if any)
p_frac <- ggplot(impact_df, aes(frac_removed)) +
  geom_histogram(bins = 40) +
  labs(x = "fraction of reads removed per sample", y = "count",
       title = sprintf("per-sample impact, thr=%s", thr_tag)) +
  theme_minimal()
ggsave(file.path(output_dir, sprintf("per_sample_fraction_removed_thr_%s.png", thr_tag)),
       p_frac, width = 5, height = 3.5, dpi = 300)

# 5) p-value histogram (only if present in current decontam version)
# p-value histogram (now always available via summary_df$p_value)
if (any(!is.na(summary_df$p_value))) {
  p_hist <- ggplot(data.frame(p = summary_df$p_value), aes(p)) +
    geom_histogram(bins = 50) +
    labs(x = p_label, y = "count", title = "decontam (prevalence) p-values") +
    theme_minimal()
  ggsave(file.path(output_dir, "decontam_score_hist.png"), p_hist, width = 5, height = 3.5, dpi = 300)
}

# 6) top genera among flagged features (useful for supplement)
if (!is.null(tax_table(ps))) {
  tx <- as.data.frame(tax_table(ps))[summary_df$feature_id, , drop = FALSE]
  tx$flagged_union <- summary_df$flagged_union
  top_gen <- sort(table(ifelse(is.na(tx$Genus), "Unclassified", tx$Genus)[tx$flagged_union]),
                  decreasing = TRUE)
  top_gen_df <- data.frame(Genus = names(top_gen), n_flagged = as.integer(top_gen))
  write.csv(top_gen_df, file = file.path(output_dir, sprintf("flagged_top_genera_thr_%s.csv", thr_tag)),
            row.names = FALSE)
}

# =============================================================================
# decision log (plain text you can use for methods/supplement)
# =============================================================================

log_lines <- c(
  "decontam prevalence filtering decision log",
  paste0("date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("feature table: ", feature_table_qza),
  paste0("taxonomy:      ", taxonomy_qza),
  paste0("metadata:      ", metadata_path),
  "",
  paste0("n negatives used: ", n_neg),
  paste0("n real samples used: ", n_smp),
  paste0("mock controls excluded from fit: ", sum(sample_data(ps)$control_role == "pos")),
  paste0("zero-read samples dropped from fit: ", length(dropped_ids)),
  "",
  "threshold sweep (see threshold_sweep.csv for full table):",
  paste(capture.output(print(sweep_res)), collapse = "\n"),
  "",
  paste0("chosen prevalence threshold: ", use_threshold, " (file tag ", thr_tag, ")"),
  paste0("union rule adds neg-only features present in >= ", min_neg_presence, " negative(s) and 0 real samples"),
  paste0("final n features removed (union): ", sum(summary_df$flagged_union)),
  paste0("% reads removed (union): ",
         round(100 * (1 - sum(sample_sums(ps_clean)) / sum(sample_sums(ps_fit))), 4)),
  paste0("n flagged features that occur in real samples: ",
         sum(summary_df$flagged_union & summary_df$sample_prev > 0)),
  paste0("n flagged features that are neg-only: ",
         sum(summary_df$flagged_union & summary_df$sample_prev == 0)),
  "",
  "notes:",
  "- prevalence bands in the scatter reflect discrete denominators (k/14 in negatives; j/80 in samples).",
  "- samples are not dropped unless they lose a very large fraction of reads post-filter (e.g., >50%).",
  "- all figures and csvs saved in this directory."
)
writeLines(log_lines, con = file.path(output_dir, paste0("DECISION_LOG_thr_", thr_tag, ".txt")))

message(sprintf("done. results written to: %s", normalizePath(output_dir)))




