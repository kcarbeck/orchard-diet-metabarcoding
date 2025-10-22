#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# Edit per-plate metadata filenames and associated tags if needed
files <- c("meta_2024A.tsv","meta_2024B.tsv","meta_2025A.tsv","meta_2025B.tsv")
tags  <- c("2024A","2024B","2025A","2025B")
stopifnot(length(files) == length(tags))

read_one <- function(path, plate_tag) {
  if (!file.exists(path)) stop(paste("Missing metadata:", path))
  md <- read_tsv(path, show_col_types = FALSE, comment = "")
  if (!"SampleID" %in% names(md) && "#SampleID" %in% names(md)) {
    names(md)[names(md) == "#SampleID"] <- "SampleID"
  }
  if (!"SampleID" %in% names(md)) stop(paste("metadata missing SampleID in", path))

  md <- md %>%
    mutate(SampleID = str_replace_all(trimws(as.character(SampleID)), "\\s+", "_")) %>%
    filter(!is.na(SampleID) & SampleID != "") %>%
    mutate(SampleID_original = SampleID, Plate = plate_tag)

  if (!"Species" %in% names(md)) stop(paste("metadata missing 'Species' column in", path))
  md
}

md_list <- Map(read_one, files, tags)
md_all  <- bind_rows(md_list)

normlab <- function(x) toupper(gsub("[^A-Za-z0-9]+","", trimws(as.character(x))))
lab <- normlab(md_all$Species)
neg <- c("EBLANK","PBLANK","BLANK","EMPTY"); pos <- "POS"
md_all$control_role <- ifelse(lab %in% pos, "pos", ifelse(lab %in% neg, "neg", "sample"))

dup_flag <- duplicated(md_all$SampleID) | duplicated(md_all$SampleID, fromLast = TRUE)
n_dup <- sum(dup_flag)

if (n_dup > 0) {
  message("WARNING: found ", n_dup, " rows with duplicated SampleID across plates.")
  md_all <- md_all %>%
    mutate(SampleID_new = paste0(Plate, "_", SampleID), SampleID = SampleID_new)
  dir.create("rename_maps", showWarnings = FALSE)
  md_all %>%
    distinct(Plate, SampleID_original, SampleID) %>%
    group_by(Plate) %>%
    arrange(SampleID_original, .by_group = TRUE) %>%
    select(SampleID = SampleID_original, new_id = SampleID) %>%
    group_split() %>%
    Map(function(df, tag) write_tsv(df, file.path("rename_maps", paste0("rename_map_", tag, ".tsv"))),
        ., unique(md_all$Plate))
} else {
  message("no duplicated SampleID detected; keeping original IDs.")
}

meta_out <- md_all %>%
  select(SampleID, control_role, Plate, everything(), -SampleID_original, -SampleID_new)

write_tsv(meta_out, "all_plates_metadata.tsv")
message("wrote: all_plates_metadata.tsv")


