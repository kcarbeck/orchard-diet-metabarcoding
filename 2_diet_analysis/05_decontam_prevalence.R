# =============================================================================
# decontam prevalence method (consolidated script)
# =============================================================================
# author: katherine carbeck
# date: 28 nov 2025
#
# this script identifies and removes contaminant sequences using the prevalence method from the decontam package
#
# how it works:
#   - contaminants are more common in negative controls than real samples
#   - decontam compares prevalence (presence/absence) between negatives and samples
#   - sequences that are more prevalent in negatives get flagged
#
# reference: davis et al. 2018 (doi:10.1186/s40168-018-0605-2)
#
# what you need:
#   1. merged feature table (.qza) from qiime2 (table_merged.qza)
#   2. merged taxonomy (.qza) from classifier (taxonomy_merged.qza)
#   3. merged metadata (.tsv) with (all_plates_metadata.tsv):
#      - SampleID column
#      - Species column (with control labels: EBLANK, PBLANK, BLANK, EMPTY, POS)
#      - Year_Plate column for batch mode (e.g., "2024_Plate01")
#
# what this script outputs:
#   - contaminant_feature_ids_thr_*.txt  (list of contaminant IDs for qiime2)
#   - contaminants_summary_thr_*.csv     (detailed summary of all features)
#   - threshold_sweep.csv                (helps you choose the right threshold)
#   - diagnostic plots                   (to check everything looks reasonable)
#   - DECISION_LOG_thr_*.txt             (for your methods section)
#
# =============================================================================

# -----------------------------------------------------------------------------
# install packages (only need to do this once)
# -----------------------------------------------------------------------------
# uncomment and run these lines if you haven't installed the packages yet:
#
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("phyloseq", "decontam"), update = FALSE, ask = FALSE)
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("jbisanz/qiime2R")

# load required packages
library(qiime2R)
library(phyloseq)
library(decontam)
library(ggplot2)

# =============================================================================
# configuration - edit these paths for your data
# =============================================================================

# input files (merged across all plates/years)
# these come from the 05_merge_plates_and_years.sh script
feature_table_qza <- "/lustre2/home/lc736_0001/orchards/merged/table_merged.qza"
taxonomy_qza      <- "/lustre2/home/lc736_0001/orchards/merged/taxonomy_merged.qza"
metadata_path     <- "/lustre2/home/lc736_0001/orchards/merged/all_plates_metadata.tsv"

# where to save output files
output_dir <- "/lustre2/home/lc736_0001/orchards/merged/decontam_output"

# -----------------------------------------------------------------------------
# control labels
# -----------------------------------------------------------------------------
# these are the values in your "Species" column that indicate controls
# case doesn't matter - we normalize everything to uppercase

# negative controls for decontam (LAB blanks only)
# these should have no real DNA, only lab/reagent contamination
# EBLANK = extraction blank, PBLANK = PCR blank, EMPTY = empty wells
negative_labels <- c("EBLANK", "PBLANK", "EMPTY")

# field blanks (substrate swabs, bag swabs, etc.)
# these capture environmental DNA from where birds forage - NOT used for decontam
# because environmental DNA is real signal, not contamination
field_blank_labels <- c("BLANK")

# positive controls (mock communities)
# these are excluded from the decontam analysis
positive_labels <- c("POS")

# -----------------------------------------------------------------------------
# column matching patterns
# -----------------------------------------------------------------------------
# regex patterns to find the right columns in your metadata
# you probably don't need to change these unless your column names are unusual

species_col_regex <- "^Species"      # column with control type info
batch_col_regex   <- "^Year_Plate$"  # unique batch identifier

# -----------------------------------------------------------------------------
# decontam parameters
# -----------------------------------------------------------------------------

# prevalence thresholds to test
# decontam will sweep through these to show you what happens at each level
# lower threshold = more conservative (flags fewer things)
# higher threshold = more aggressive (flags more things)
thr_grid <- c(0.05, 0.10, 0.20, 0.30, 0.50)

# neg-only rule: flag features that appear in negatives but never in real samples
# these are almost certainly contaminants
# setting this to 1 means: flag if in >= 1 negative and 0 real samples
min_neg_presence <- 1

# =============================================================================
# helper functions
# =============================================================================

# normalize labels to uppercase, remove special characters
# this makes matching more robust (e.g., "E-Blank" matches "EBLANK")
normlab <- function(x) {
  toupper(gsub("[^A-Za-z0-9]+", "", trimws(as.character(x))))
}

# convert threshold to filename-friendly string (0.10 -> "010")
thr_tag_fn <- function(x) {
  gsub("\\.", "", sprintf("%.2f", x))
}

# fix newlines inside quoted cells (common problem with excel exports)
# this prevents weird parsing errors
fix_newlines_inside_quotes <- function(path) {
  txt <- readChar(path, file.info(path)$size, useBytes = TRUE)
  txt <- gsub("\r\n?", "\n", txt, perl = TRUE)
  ch <- strsplit(txt, "", fixed = TRUE)[[1]]
  out <- ch
  inq <- FALSE
  for (i in seq_along(ch)) {
    if (ch[i] == "\"") inq <- !inq
    if (ch[i] == "\n" && inq) out[i] <- " "
  }
  tmp <- tempfile(fileext = ".tsv")
  writeChar(paste(out, collapse = ""), tmp, eos = NULL, useBytes = TRUE)
  tmp
}

# prepare metadata file for qza_to_phyloseq
prep_metadata <- function(path) {
  fixed_path <- fix_newlines_inside_quotes(path)
  
  # try qiime2R first, fall back to standard read if no #q2:types line
  tryCatch({
    md <- qiime2R::read_q2metadata(fixed_path)
  }, error = function(e) {
    if (grepl("q2:types", e$message, ignore.case = TRUE)) {
      message("  metadata doesn't have #q2:types line, reading as standard TSV...")
      md <<- read.delim(fixed_path, sep = "\t", header = TRUE, 
                        stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      stop(e)
    }
  })
  
  # handle different column name formats
  if (!("SampleID" %in% names(md)) && "#SampleID" %in% names(md)) {
    names(md)[names(md) == "#SampleID"] <- "SampleID"
  }
  if (!("SampleID" %in% names(md))) {
    stop("metadata must include SampleID column")
  }
  
  # clean up sample ids (remove whitespace, replace spaces with underscores)
  md$SampleID <- gsub("\\s+", "_", trimws(as.character(md$SampleID)))
  md <- md[!is.na(md$SampleID) & md$SampleID != "", , drop = FALSE]
  
  # warn about duplicates (keep first occurrence)
  if (anyDuplicated(md$SampleID)) {
    dups <- unique(md$SampleID[duplicated(md$SampleID)])
    message("warning: duplicate sample ids removed (kept first): ", 
            paste(dups, collapse = ", "))
    md <- md[!duplicated(md$SampleID), , drop = FALSE]
  }
  
  # write cleaned metadata to temp file
  tmp <- tempfile(fileext = ".tsv")
  write.table(md, tmp, sep = "\t", quote = FALSE, row.names = FALSE)
  tmp
}

# =============================================================================
# load data and build phyloseq object
# =============================================================================

# create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("loading data...")

# prepare metadata and load into phyloseq
meta_clean <- prep_metadata(metadata_path)
ps <- qza_to_phyloseq(
  features = feature_table_qza, 
  taxonomy = taxonomy_qza, 
  metadata = meta_clean
)

message("  loaded ", nsamples(ps), " samples and ", ntaxa(ps), " features")
#   loaded 448 samples and 40294 features

# -----------------------------------------------------------------------------
# identify controls
# -----------------------------------------------------------------------------
sd <- sample_data(ps)

# find the species and batch columns
species_col <- grep(species_col_regex, names(sd), ignore.case = TRUE, value = TRUE)[1]
batch_col <- grep(batch_col_regex, names(sd), ignore.case = TRUE, value = TRUE)[1]

if (is.na(species_col)) {
  stop("could not find Species column in metadata (adjust species_col_regex)")
}
if (is.na(batch_col)) {
  stop("could not find Year_Plate column in metadata (adjust batch_col_regex)")
}

# assign control roles based on species column
labs <- normlab(sd[[species_col]])
sd$control_role <- factor(
  ifelse(labs %in% normlab(positive_labels), "pos",
         ifelse(labs %in% normlab(negative_labels), "neg",
                ifelse(labs %in% normlab(field_blank_labels), "field", "sample"))),
  levels = c("sample", "neg", "field", "pos")
)
sd$is_neg <- sd$control_role == "neg"
sd$plate_batch <- factor(as.character(sd[[batch_col]]))
sample_data(ps) <- sd

# print control counts
message("\ncontrol counts:")
print(table(sd$control_role, useNA = "ifany"))
message("  note: 'neg' = lab blanks (used for decontam)")
message("        'field' = field blanks (excluded from decontam - environmental DNA is real signal)")

# save control role mapping for reference
write.csv(
  data.frame(SampleID = rownames(sd), control_role = as.character(sd$control_role)),
  file.path(output_dir, "sample_control_roles.csv"), 
  row.names = FALSE
)
#control counts:
#sample    neg  field    pos
#   403     35      2      8
#  note: 'neg' = lab blanks (used for decontam)
#        'field' = field blanks (excluded from decontam - environmental DNA is real signal)


# -----------------------------------------------------------------------------
# prepare for decontam
# -----------------------------------------------------------------------------
# exclude positive controls (mocks) AND field blanks from the analysis
# decontam is designed for lab negatives vs real samples only
# field blanks contain environmental DNA which is real signal, not contamination
ps_fit <- prune_samples(sd$control_role %in% c("sample", "neg"), ps)

# remove samples with zero reads (decontam will ignore them anyway)
zs <- sample_sums(ps_fit)
dropped_ids <- names(zs)[zs == 0]
if (length(dropped_ids)) {
  message("\ndropped ", length(dropped_ids), " zero-read samples")
  write.table(
    data.frame(SampleID = dropped_ids),
    file.path(output_dir, "dropped_zero_read_samples.txt"),
    quote = FALSE, sep = "\t", row.names = FALSE
  )
}
ps_fit <- prune_samples(sample_sums(ps_fit) > 0, ps_fit)

# get counts for prevalence calculations
counts <- as(otu_table(ps_fit), "matrix")
if (!taxa_are_rows(ps_fit)) counts <- t(counts)
counts_bin <- counts > 0  # convert to presence/absence

neg_flag <- sample_data(ps_fit)$is_neg
n_neg <- sum(neg_flag)
n_smp <- sum(!neg_flag)

message("\nusing ", n_neg, " negative controls and ", n_smp, " real samples")

# find features that only appear in negatives (strong evidence of contamination)
neg_prev_counts <- rowSums(counts_bin[, neg_flag, drop = FALSE])
smp_prev_counts <- rowSums(counts_bin[, !neg_flag, drop = FALSE])
neg_only_ids <- names(which((neg_prev_counts >= min_neg_presence) & (smp_prev_counts == 0)))

message("  ", length(neg_only_ids), " features appear only in negatives (will be flagged)")

# dropped 18 zero-read samples
# using 23 negative controls and 397 real samples
#  38 features appear only in negatives (will be flagged)


# =============================================================================
# threshold sweep
# =============================================================================
# test multiple thresholds to see how many features/reads get removed
# this helps you choose an appropriate threshold

message("\nrunning threshold sweep...")

sweep_fn <- function(t) {
  # run decontam at this threshold
  ct <- isContaminant(ps_fit, method = "prevalence", neg = "is_neg", 
                      threshold = t, batch = "plate_batch")
  
  # combine with neg-only rule
  ids_union <- union(rownames(ct)[ct$contaminant], neg_only_ids)
  
  # calculate impact
  ps_tmp <- prune_taxa(!(taxa_names(ps_fit) %in% ids_union), ps_fit)
  
  data.frame(
    threshold = t,
    n_features_removed = length(ids_union),
    pct_reads_removed = 100 * (1 - sum(sample_sums(ps_tmp)) / sum(sample_sums(ps_fit)))
  )
}

sweep_res <- do.call(rbind, lapply(thr_grid, sweep_fn))
write.csv(sweep_res, file.path(output_dir, "threshold_sweep.csv"), row.names = FALSE)

message("\nthreshold sweep results:")
print(sweep_res)

# plot the sweep results
p_sweep <- ggplot(sweep_res, aes(threshold, pct_reads_removed)) +
  geom_line(linewidth = 1) + 
  geom_point(size = 3) +
  labs(
    x = "prevalence threshold",
    y = "% reads removed",
    title = "decontam threshold sweep"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "threshold_sweep.png"), p_sweep, 
       width = 5, height = 3.5, dpi = 300)

# threshold sweep results:
#  threshold n_features_removed pct_reads_removed
# 1      0.05                 49         0.2467184
# 2      0.10                 76         0.7520995
# 3      0.20                100         2.1765559
# 4      0.30                112         7.1294023
# 5      0.50                126        12.2777405
# looks like between 0.10 and 0.20 is a good threshold, but want to check what those asvs are (next section)

# =============================================================================
# investigate borderline features 
# =============================================================================
# this section helps you decide between thresholds by showing what features would be additionally removed at higher thresholds

# compare two thresholds (0.10 vs 0.20)
thr_low <- 0.10
thr_high <- 0.20

ct_low <- isContaminant(ps_fit, method = "prevalence", neg = "is_neg",
                        threshold = thr_low, batch = "plate_batch")
ct_high <- isContaminant(ps_fit, method = "prevalence", neg = "is_neg",
                         threshold = thr_high, batch = "plate_batch")

# features flagged at higher threshold but NOT at lower threshold
ids_low <- union(rownames(ct_low)[ct_low$contaminant], neg_only_ids)
ids_high <- union(rownames(ct_high)[ct_high$contaminant], neg_only_ids)
borderline_ids <- setdiff(ids_high, ids_low)

message("\n", length(borderline_ids), " features flagged at ", thr_high, 
        " but not at ", thr_low, ":")
# 24 features flagged at 0.2 but not at 0.1:

# get taxonomy for borderline features
tax <- as.data.frame(tax_table(ps_fit))
borderline_tax <- tax[borderline_ids, , drop = FALSE]

# calculate prevalence stats (proportion of samples each feature appears in)
neg_prev <- rowMeans(counts_bin[, neg_flag, drop = FALSE])
smp_prev <- rowMeans(counts_bin[, !neg_flag, drop = FALSE])

# get prevalence stats for borderline features
borderline_df <- data.frame(
  feature_id = borderline_ids,
  neg_prev = neg_prev[borderline_ids],
  sample_prev = smp_prev[borderline_ids],
  neg_reads = rowSums(counts[borderline_ids, neg_flag, drop = FALSE]),
  sample_reads = rowSums(counts[borderline_ids, !neg_flag, drop = FALSE]),
  borderline_tax[borderline_ids, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
)

# sort by sample reads (most impactful first)
borderline_df <- borderline_df[order(-borderline_df$sample_reads), ]

# print top features you'd be removing
message("\ntop 20 borderline features by sample reads:")
print(head(borderline_df, 20))

# save full list for review
write.csv(borderline_df, 
          file.path(output_dir, paste0("borderline_features_", thr_low, "_vs_", thr_high, ".csv")),
          row.names = FALSE)
message("full list saved to: borderline_features_", thr_low, "_vs_", thr_high, ".csv")

#took a look at the borderline features and they don't look like contaminants-they're all reasonable prey species, so we'll use the 0.10 threshold

# =============================================================================
# apply chosen threshold
# =============================================================================
use_threshold <- 0.10

message("\napplying threshold: ", use_threshold)
thr_tag <- thr_tag_fn(use_threshold)

# run decontam with chosen threshold
ct_use <- isContaminant(ps_fit, method = "prevalence", neg = "is_neg",
                        threshold = use_threshold, batch = "plate_batch")

# combine decontam results with neg-only rule
contam_ids <- union(rownames(ct_use)[ct_use$contaminant], neg_only_ids)

message("  flagged ", length(contam_ids), " features as contaminants")
# applying threshold: 0.1
#  flagged 76 features as contaminants

# -----------------------------------------------------------------------------
# create summary table
# -----------------------------------------------------------------------------
neg_prev <- rowMeans(counts_bin[, neg_flag, drop = FALSE])
smp_prev <- rowMeans(counts_bin[, !neg_flag, drop = FALSE])

summary_df <- data.frame(
  feature_id = rownames(counts),
  neg_prev = neg_prev[rownames(counts)],
  sample_prev = smp_prev[rownames(counts)],
  neg_reads = rowSums(counts[, neg_flag, drop = FALSE])[rownames(counts)],
  sample_reads = rowSums(counts[, !neg_flag, drop = FALSE])[rownames(counts)],
  contaminant = rownames(counts) %in% contam_ids,
  stringsAsFactors = FALSE
)

write.csv(summary_df, 
          file.path(output_dir, paste0("contaminants_summary_thr_", thr_tag, ".csv")),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# create detailed table of flagged contaminants with taxonomy
# -----------------------------------------------------------------------------
tax_all <- as.data.frame(tax_table(ps_fit))
contam_tax <- tax_all[contam_ids, , drop = FALSE]

contam_details <- data.frame(
  feature_id = contam_ids,
  neg_prev = neg_prev[contam_ids],
  sample_prev = smp_prev[contam_ids],
  neg_reads = rowSums(counts[contam_ids, neg_flag, drop = FALSE]),
  sample_reads = rowSums(counts[contam_ids, !neg_flag, drop = FALSE]),
  contam_tax[contam_ids, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")],
  stringsAsFactors = FALSE
)

# sort by sample_reads (most impactful first)
contam_details <- contam_details[order(-contam_details$sample_reads), ]

write.csv(contam_details, 
          file.path(output_dir, paste0("contaminants_with_taxonomy_thr_", thr_tag, ".csv")),
          row.names = FALSE)

message("  saved contaminant details with taxonomy: contaminants_with_taxonomy_thr_", thr_tag, ".csv")


# -----------------------------------------------------------------------------
# save contaminant id list for qiime2
# -----------------------------------------------------------------------------
# this file is used by qiime feature-table filter-features
write.table(
  data.frame(FeatureID = contam_ids),
  file.path(output_dir, paste0("contaminant_feature_ids_thr_", thr_tag, ".txt")),
  quote = FALSE, sep = "\t", row.names = FALSE
)

# =============================================================================
# field blank analysis
# =============================================================================
# field blanks capture environmental DNA from sampling substrates
# this section summarizes what's in them and flags those taxa in the output

message("\nanalyzing field blanks...")

# get field blank samples from the original phyloseq object (before we excluded them)
field_samples <- rownames(sample_data(ps))[sample_data(ps)$control_role == "field"]
n_field <- length(field_samples)

if (n_field > 0) {
  message("  found ", n_field, " field blank samples")
  
  # subset to field blanks only
  ps_field <- prune_samples(sample_names(ps) %in% field_samples, ps)
  ps_field <- prune_taxa(taxa_sums(ps_field) > 0, ps_field)  # remove absent taxa
  
  # get counts matrix for field blanks
  field_counts <- as(otu_table(ps_field), "matrix")
  if (!taxa_are_rows(ps_field)) field_counts <- t(field_counts)
  
  # calculate field blank stats
  field_taxa <- rownames(field_counts)
  field_prev <- rowSums(field_counts > 0) / n_field  # prevalence in field blanks
  field_reads <- rowSums(field_counts)  # total reads in field blanks
  
  # get taxonomy
  tax_field <- as.data.frame(tax_table(ps_field))
  
  # create field blank summary
  field_summary <- data.frame(
    feature_id = field_taxa,
    field_prev = field_prev[field_taxa],
    field_reads = field_reads[field_taxa],
    tax_field[field_taxa, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")],
    stringsAsFactors = FALSE
  )
  
  # sort by reads (most abundant first)
  field_summary <- field_summary[order(-field_summary$field_reads), ]
  
  write.csv(field_summary,
            file.path(output_dir, "field_blank_taxa_summary.csv"),
            row.names = FALSE)
  
  message("  ", length(field_taxa), " taxa detected in field blanks")
  message("  saved: field_blank_taxa_summary.csv")
  
  # print top taxa in field blanks
  message("\n  top 10 taxa in field blanks by read count:")
  print(head(field_summary[, c("field_reads", "Order", "Family", "Genus", "Species")], 10))
  
  # -----------------------------------------------------------------------------
  # add field blank flag to main summary table
  # -----------------------------------------------------------------------------
  # reload the summary and add field blank column
  summary_df$in_field_blank <- summary_df$feature_id %in% field_taxa
  
  # also calculate how many field blank reads each feature has
  summary_df$field_blank_reads <- 0
  summary_df$field_blank_reads[match(field_taxa, summary_df$feature_id)] <- field_reads[field_taxa]
  
  # re-save the updated summary
  write.csv(summary_df, 
            file.path(output_dir, paste0("contaminants_summary_thr_", thr_tag, ".csv")),
            row.names = FALSE)
  
  message("\n  updated contaminants_summary with 'in_field_blank' and 'field_blank_reads' columns")
  message("  ", sum(summary_df$in_field_blank), " features appear in field blanks")
  message("  ", sum(summary_df$in_field_blank & !summary_df$contaminant), 
          " of these are NOT flagged as contaminants (kept in final data)")
  
} else {
  message("  no field blanks found - skipping field blank analysis")
}

# updated contaminants_summary with 'in_field_blank' and 'field_blank_reads' columns
#  5 features appear in field blanks
#  5 of these are NOT flagged as contaminants (kept in final data)

# =============================================================================
# diagnostic plots
# =============================================================================

# -----------------------------------------------------------------------------
# plot 1: library sizes by control role
# -----------------------------------------------------------------------------
# negatives should have much lower read counts than real samples
# if negatives have high counts, something may be wrong

lib_df <- data.frame(
  sample_id = sample_names(ps_fit),
  control_role = as.character(sample_data(ps_fit)$control_role),
  reads = as.numeric(sample_sums(ps_fit))
)
pal <- c(sample = "#0072B2", neg = "#D55E00", pos = "#009E73")

p_lib <- ggplot(lib_df, aes(control_role, reads, fill = control_role)) +
  geom_violin(alpha = 0.8, trim = FALSE) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "control role",
    y = "reads",
    title = "library sizes by control role"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "library_sizes_by_role.png"), p_lib, 
       width = 5, height = 4, dpi = 300)

# -----------------------------------------------------------------------------
# plot 2: prevalence scatter
# -----------------------------------------------------------------------------
# this shows prevalence in samples vs negatives for each feature
# contaminants should be in the upper left (high in negatives, low in samples)

plot_df <- data.frame(
  sample_prev = smp_prev,
  neg_prev = neg_prev,
  contaminant = names(smp_prev) %in% contam_ids
)

p_prev <- ggplot(plot_df, aes(sample_prev, neg_prev, color = contaminant)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_manual(
    values = c("FALSE" = "#0072B2", "TRUE" = "#D55E00"),
    labels = c("kept", "flagged"),
    name = "status"
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "sample prevalence",
    y = "negative prevalence",
    title = paste0("prevalence scatter (threshold=", use_threshold, ")")
  ) +
  theme_minimal()

ggsave(file.path(output_dir, paste0("prevalence_scatter_thr_", thr_tag, ".png")),
       p_prev, width = 5, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# plot 3: per-sample impact
# -----------------------------------------------------------------------------
# shows what fraction of reads were removed from each sample
# most samples should have very low fractions removed

ps_clean <- prune_taxa(!(taxa_names(ps_fit) %in% contam_ids), ps_fit)
impact_df <- data.frame(
  sample_id = sample_names(ps_fit),
  reads_before = as.integer(sample_sums(ps_fit)),
  reads_after = as.integer(sample_sums(ps_clean)),
  frac_removed = 1 - sample_sums(ps_clean) / pmax(1, sample_sums(ps_fit))
)

p_impact <- ggplot(impact_df, aes(frac_removed)) +
  geom_histogram(bins = 40, fill = "#0072B2", alpha = 0.8) +
  labs(
    x = "fraction of reads removed",
    y = "count",
    title = paste0("per-sample impact (threshold=", use_threshold, ")")
  ) +
  theme_minimal()

ggsave(file.path(output_dir, paste0("per_sample_impact_thr_", thr_tag, ".png")),
       p_impact, width = 5, height = 3.5, dpi = 300)

# =============================================================================
# decision log
# =============================================================================
# this creates a text file you can use for your methods section

log_lines <- c(
  "DECONTAM PREVALENCE FILTERING - DECISION LOG",
  paste0("date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "=== input files ===",
  paste0("feature table: ", feature_table_qza),
  paste0("taxonomy: ", taxonomy_qza),
  paste0("metadata: ", metadata_path),
  "",
  "=== sample counts ===",
  paste0("negative controls: ", n_neg),
  paste0("real samples: ", n_smp),
  paste0("batches (Year_Plate): ", paste(levels(sample_data(ps_fit)$plate_batch), collapse = ", ")),
  paste0("zero-read samples dropped: ", length(dropped_ids)),
  "",
  "=== threshold sweep ===",
  paste(capture.output(print(sweep_res)), collapse = "\n"),
  "",
  "=== chosen parameters ===",
  paste0("prevalence threshold: ", use_threshold),
  paste0("neg-only rule: features in >= ", min_neg_presence, " negatives and 0 real samples"),
  "",
  "=== results ===",
  paste0("total features flagged: ", length(contam_ids)),
  paste0("% reads removed: ", round(100 * (1 - sum(sample_sums(ps_clean)) / sum(sample_sums(ps_fit))), 4)),
  "",
  "=== reference ===",
  "davis et al. 2018. simple statistical identification and removal of contaminant",
  "sequences in marker-gene and metagenomics data. microbiome 6:226.",
  "doi:10.1186/s40168-018-0605-2"
)

writeLines(log_lines, file.path(output_dir, paste0("DECISION_LOG_thr_", thr_tag, ".txt")))

# =============================================================================
# save r objects for later use
# =============================================================================
saveRDS(
  list(ps_fit = ps_fit, ps_clean = ps_clean, contam_ids = contam_ids, sweep = sweep_res),
  file.path(output_dir, paste0("decontam_objects_thr_", thr_tag, ".rds"))
)

# save session info for reproducibility
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))


# ============================================
#              DECONTAM DONE :) 
# results saved to: /lustre2/home/lc736_0001/orchards/merged/decontam_output
# ============================================
# 
# next steps:
#   1. review threshold_sweep.csv and diagnostic plots
#   2. if threshold looks good, use contaminant_feature_ids_thr_010.txt
#   3. filter your table with qiime feature-table filter-features (06_qc_and_get_analysis_files.sh)
