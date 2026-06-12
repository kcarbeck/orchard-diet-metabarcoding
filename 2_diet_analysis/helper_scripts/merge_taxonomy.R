#!/usr/bin/env Rscript
# =============================================================================
# merge taxonomy files from multiple plates/years
# =============================================================================
# author: katherine carbeck
# date: nov 2025
#
# this is a HELPER SCRIPT called by 04.2_merge_plates_and_years.sh
# it merges taxonomy from exported qiime2 files
#
# when the same asv appears in multiple plates, we keep:
#   1. the classification with the deepest rank (most specific)
#   2. on ties, the one with highest confidence score
#
# input: exported taxonomy tsv files from qiime tools export
# output: taxonomy_merged.tsv (ready to import back to qiime2)
#
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# =============================================================================
# configuration - edit these paths if needed
# =============================================================================

# list of exported taxonomy tsv files
# these are created by: qiime tools export --input-path taxonomy.qza --output-path export_taxonomy/
# the actual file will be at export_taxonomy/taxonomy.tsv
files <- c(
  "export_nb_classified_taxonomy_111125/taxonomy.tsv",
  "export_nb_classified_taxonomy_141025/taxonomy.tsv",
  "export_nb_classified_taxonomy_112425/taxonomy.tsv"
)

# output filename
output_file <- "taxonomy_merged.tsv"

# =============================================================================
# merge taxonomy
# =============================================================================

# keep only files that exist
files <- files[file.exists(files)]

if (length(files) == 0) {
  stop("no taxonomy tsv files found - check file paths in configuration")
}

message("found ", length(files), " taxonomy file(s) to merge:")
for (f in files) {
  message("  - ", f)
}

# read and combine all taxonomy files
tx <- bind_rows(lapply(files, read_tsv, show_col_types = FALSE))

# check for required columns
required_cols <- c("Feature ID", "Taxon", "Confidence")
if (!all(required_cols %in% names(tx))) {
  missing <- setdiff(required_cols, names(tx))
  stop("taxonomy files missing required columns: ", paste(missing, collapse = ", "))
}

message("\ntotal rows before merging: ", nrow(tx))
message("unique features: ", length(unique(tx$`Feature ID`)))

# -----------------------------------------------------------------------------
# calculate taxonomy depth
# -----------------------------------------------------------------------------
# depth = number of taxonomic ranks assigned
# e.g., "k__Animalia;p__Arthropoda;c__Insecta" has depth 3

depth <- function(t) {
  ifelse(is.na(t) | t == "", 0L, str_count(t, ";") + 1L)
}

# -----------------------------------------------------------------------------
# merge: keep deepest classification, then highest confidence
# -----------------------------------------------------------------------------
tx_merged <- tx %>%
  group_by(`Feature ID`) %>%
  arrange(desc(depth(Taxon)), desc(Confidence)) %>%
  slice(1) %>%
  ungroup()

message("rows after merging: ", nrow(tx_merged))

# write output
write_tsv(tx_merged[, c("Feature ID", "Taxon", "Confidence")], output_file)
message("\nwrote: ", output_file)