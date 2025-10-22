#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# Edit list if needed; these are the exported per-plate taxonomy TSVs
files <- c(
  "export_taxonomy_2024A/taxonomy.tsv",
  "export_taxonomy_2024B/taxonomy.tsv",
  "export_taxonomy_2025A/taxonomy.tsv",
  "export_taxonomy_2025B/taxonomy.tsv"
)

files <- files[file.exists(files)]
if (length(files) == 0) {
  stop("No taxonomy TSVs found. Expected export_taxonomy_*/taxonomy.tsv")
}

tx <- bind_rows(lapply(files, read_tsv, show_col_types = FALSE))

if (!all(c("Feature.ID","Taxon","Confidence") %in% names(tx))) {
  stop("taxonomy.tsv must contain Feature.ID, Taxon, Confidence columns")
}

depth <- function(t) ifelse(is.na(t) | t == "", 0L, str_count(t, ";") + 1L)

tx_merged <- tx %>%
  group_by(Feature.ID) %>%
  arrange(desc(depth(Taxon)), desc(Confidence)) %>%
  slice(1) %>%
  ungroup()

write_tsv(tx_merged[, c("Feature.ID","Taxon","Confidence")], "taxonomy_merged.tsv")
message("wrote taxonomy_merged.tsv (n=", nrow(tx_merged), ")")


