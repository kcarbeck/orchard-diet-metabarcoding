# =============================================================================
# merge metadata across plates and years
# =============================================================================
# author: katherine carbeck
# date: nov 2025
#
# this script combines metadata from multiple sequencing plates/years into
# one unified file for downstream analysis
#
# what this script does:
#   1. reads metadata from each plate/year
#   2. normalizes column names (handles #SampleID, sample-id, etc.)
#   3. creates a Year_Plate column for batch identification in decontam
#   4. detects duplicate sample ids (usually just controls)
#   5. prefixes duplicate control ids to make them unique
#   6. generates rename maps for qiime2 if needed
#
# output files:
#   - all_plates_metadata.tsv        merged metadata with unique ids
#   - rename_maps/*.tsv              rename mappings for qiime2 (if duplicates found)
#   - duplicate_controls_log.txt     log of which ids were renamed
#
# =============================================================================

# load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# =============================================================================
# configuration - edit these for your data
# =============================================================================

# list of input metadata files
# format: list of list(path = "filename", year = "YYYY")
# if a file contains multiple plates, they will be extracted from PlateName column
input_files <- list(
  list(path = "/lustre2/home/lc736_0001/orchards/merged/input/2024/2024_metadata_3.txt", year = "2024"),
  list(path = "/lustre2/home/lc736_0001/orchards/merged/input/2025_p1/2025_fecal_sampling_metadata_plate1.txt", year = "2025"),
  list(path = "/lustre2/home/lc736_0001/orchards/merged/input/2025_p2/2025_fecal_sampling_metadata_plate2.txt", year = "2025")
)

# output files
output_metadata <- "/lustre2/home/lc736_0001/orchards/merged/all_plates_metadata.tsv"
output_dir <- "rename_maps"

# control labels (case-insensitive)
# these are the values in the "Species" column that indicate controls
negative_labels <- c("EBLANK", "PBLANK", "BLANK", "EMPTY")
positive_labels <- c("POS")

# =============================================================================
# helper functions
# =============================================================================

# normalize labels to uppercase, remove special characters
# this makes matching more robust
normlab <- function(x) {
  toupper(gsub("[^A-Za-z0-9]+", "", trimws(as.character(x))))
}

# -----------------------------------------------------------------------------
# read a single metadata file
# -----------------------------------------------------------------------------
# handles different column name formats and creates Year_Plate identifier

read_metadata <- function(path, year) {
  
  # check file exists
  if (!file.exists(path)) {
    stop(paste("file not found:", path))
  }
  
  message("  reading: ", path)
  
  # read the file (all columns as character to avoid type issues)
  md <- read_tsv(path, show_col_types = FALSE, comment = "", 
                 col_types = cols(.default = "c"))
  
  # remove the #q2:types row if present (qiime2 metadata format)
  if (any(grepl("^#q2:types$", md[[1]], ignore.case = TRUE))) {
    md <- md[!grepl("^#q2:types$", md[[1]], ignore.case = TRUE), ]
  }
  
  # find and normalize the SampleID column
  # different tools use different names: #SampleID, sample-id, SampleID
  id_col <- grep("^(#?SampleID|sample-id)$", names(md), ignore.case = TRUE, value = TRUE)
  if (!length(id_col)) {
    stop(paste("could not find SampleID column in", path))
  }
  names(md)[names(md) == id_col[1]] <- "SampleID"
  
  # find the Species column (needed for identifying controls)
  species_col <- grep("^Species", names(md), ignore.case = TRUE, value = TRUE)
  if (!length(species_col)) {
    stop(paste("could not find Species column in", path))
  }
  names(md)[names(md) == species_col[1]] <- "Species"
  
  # find or create PlateName column
  plate_col <- grep("^PlateName$", names(md), ignore.case = TRUE, value = TRUE)
  if (!length(plate_col)) {
    # try to infer plate from filename (e.g., "plate1" -> "Plate01")
    plate_match <- str_extract(basename(path), "[Pp]late\\d+")
    if (!is.na(plate_match)) {
      plate_num <- str_extract(plate_match, "\\d+")
      md$PlateName <- paste0("Plate", str_pad(plate_num, 2, pad = "0"))
      message("    inferred PlateName from filename: ", md$PlateName[1])
    } else {
      stop(paste("could not determine plate from", path))
    }
  }
  
  # clean up sample ids
  # remove leading/trailing whitespace, replace internal spaces with underscores
  md <- md %>%
    mutate(
      SampleID = str_replace_all(trimws(as.character(SampleID)), "\\s+", "_"),
      SampleID_original = SampleID,  # keep original for reference
      Year = year
    ) %>%
    filter(!is.na(SampleID) & SampleID != "")  # remove empty rows
  
  # standardize PlateName format (Plate01, Plate02, etc.)
  md <- md %>%
    mutate(PlateName = case_when(
      grepl("^Plate0\\d$", PlateName) ~ PlateName,  # already formatted
      grepl("^Plate\\d$", PlateName) ~ paste0("Plate0", str_extract(PlateName, "\\d")),
      TRUE ~ PlateName
    ))
  
  # create Year_Plate batch identifier
  # this is unique across years (e.g., "2024_Plate01" vs "2025_Plate01")
  # used by decontam for batch mode
  md <- md %>%
    mutate(Year_Plate = paste0(Year, "_", PlateName))
  
  message("    found ", nrow(md), " samples in ", length(unique(md$Year_Plate)), " plate(s)")
  
  md
}

# =============================================================================
# read and merge all metadata files
# =============================================================================

message("\nreading metadata files...")
md_list <- lapply(input_files, function(x) read_metadata(x$path, x$year))
md_all <- bind_rows(md_list)

message("\ntotal rows: ", nrow(md_all))
message("plates: ", paste(unique(md_all$Year_Plate), collapse = ", "))
#reading metadata files...
#  reading: /lustre2/home/lc736_0001/orchards/merged/input/2024/2024_metadata_3.txt
#    found 256 samples in 3 plate(s)
#  reading: /lustre2/home/lc736_0001/orchards/merged/input/2025_p1/2025_fecal_sampling_metadata_plate1.txt
#    found 96 samples in 1 plate(s)
#  reading: /lustre2/home/lc736_0001/orchards/merged/input/2025_p2/2025_fecal_sampling_metadata_plate2.txt
#    found 96 samples in 1 plate(s)
#
# total rows: 448
# plates: 2024_Plate03, 2024_Plate02, 2024_Plate01, 2025_Plate01, 2025_Plate02


# =============================================================================
# assign control roles
# =============================================================================
# control_role will be used to filter out controls before analysis

lab <- normlab(md_all$Species)
md_all$control_role <- ifelse(
  lab %in% normlab(positive_labels), "pos",
  ifelse(lab %in% normlab(negative_labels), "neg", "sample")
)

message("\ncontrol role counts:")
print(table(md_all$control_role, useNA = "ifany"))
# control role counts:
# neg    pos sample
# 37      8    403


# =============================================================================
# find duplicate sample ids
# =============================================================================

# find all ids that appear more than once
dup_ids <- md_all$SampleID[duplicated(md_all$SampleID) | duplicated(md_all$SampleID, fromLast = TRUE)]
dup_ids_unique <- unique(dup_ids)

# separate duplicates into controls vs real samples
dup_control_rows <- md_all$SampleID %in% dup_ids_unique & md_all$control_role != "sample"
dup_sample_rows <- md_all$SampleID %in% dup_ids_unique & md_all$control_role == "sample"

n_dup_controls <- sum(dup_control_rows)
n_dup_samples <- sum(dup_sample_rows)

message("\nduplicate analysis:")
message("  duplicate control ids: ", n_dup_controls, " rows")
message("  duplicate real sample ids: ", n_dup_samples, " rows")

# warn if real samples have duplicates (this is usually a problem)
if (n_dup_samples > 0) {
  warning("found duplicate REAL sample ids - these should be investigated!")
  dup_samples <- md_all %>% 
    filter(SampleID %in% dup_ids_unique, control_role == "sample") %>%
    select(SampleID, Year_Plate, Species)
  print(dup_samples)
}

#duplicate analysis:
#  duplicate control ids: 18 rows
#  duplicate real sample ids: 0 rows


# =============================================================================
# prefix duplicate control ids
# =============================================================================
# only controls get prefixed - real sample ids stay as-is
# format: Year_Plate_OriginalID (e.g., "2024_Plate02_P1")

if (n_dup_controls > 0) {
  message("\nprefixing duplicate control ids with Year_Plate...")
  
  # create new ids for duplicate controls only
  md_all <- md_all %>%
    mutate(
      SampleID_new = ifelse(
        SampleID %in% dup_ids_unique & control_role != "sample",
        paste0(Year_Plate, "_", SampleID),
        SampleID
      )
    )
  
  # log the renames
  rename_log <- md_all %>%
    filter(SampleID != SampleID_new) %>%
    select(original_id = SampleID, new_id = SampleID_new, Year_Plate, Species, control_role) %>%
    arrange(Year_Plate, original_id)
  
  message("  renamed ", nrow(rename_log), " control ids")
  
  # save the log
  write_tsv(rename_log, "duplicate_controls_log.txt")
  message("  wrote: duplicate_controls_log.txt")
  
  # apply the new ids
  md_all$SampleID <- md_all$SampleID_new
  
  # -------------------------------------------------------------------------
  # generate rename maps for qiime2
  # -------------------------------------------------------------------------
  # qiime2 needs these to rename sample ids in feature tables
  # use with: qiime feature-table rename-ids
  
  dir.create(output_dir, showWarnings = FALSE)
  
  plates_with_renames <- unique(rename_log$Year_Plate)
  
  for (yp in plates_with_renames) {
    plate_renames <- rename_log %>%
      filter(Year_Plate == yp) %>%
      select(SampleID = original_id, new_id)
    
    outfile <- file.path(output_dir, paste0("rename_map_", yp, ".tsv"))
    write_tsv(plate_renames, outfile)
    message("  wrote: ", outfile, " (", nrow(plate_renames), " renames)")
  }
  
} else {
  message("\nno duplicate control ids detected - keeping original ids")
  md_all$SampleID_new <- md_all$SampleID
}
#prefixing duplicate control ids with Year_Plate...
#  renamed 18 control ids
#  wrote: duplicate_controls_log.txt
#  wrote: rename_maps/rename_map_2024_Plate02.tsv (5 renames)
#  wrote: rename_maps/rename_map_2024_Plate03.tsv (4 renames)
#  wrote: rename_maps/rename_map_2025_Plate01.tsv (8 renames)
#  wrote: rename_maps/rename_map_2025_Plate02.tsv (1 renames)


# =============================================================================
# write final output
# =============================================================================

# select and order columns for output
essential_cols <- c("SampleID", "control_role", "Year", "PlateName", "Year_Plate")
other_cols <- setdiff(names(md_all), c(essential_cols, "SampleID_original", "SampleID_new"))

meta_out <- md_all %>%
  select(all_of(essential_cols), all_of(other_cols)) %>%
  arrange(Year, PlateName, SampleID)

# write the merged metadata
write_tsv(meta_out, output_metadata)
message("\nwrote: ", output_metadata, " (", nrow(meta_out), " samples)")
# wrote: /lustre2/home/lc736_0001/orchards/merged/all_plates_metadata.tsv (448 samples)

# =============================================================================
# summary
# =============================================================================

message("\n============================================")
message("merge complete!")
message("============================================")
message("\ntotal samples: ", nrow(meta_out))
message("year_plate batches: ", length(unique(meta_out$Year_Plate)))

cat("\nsamples per Year_Plate:\n")
print(table(meta_out$Year_Plate))

cat("\ncontrol roles per Year_Plate:\n")
print(table(meta_out$Year_Plate, meta_out$control_role))
#               neg pos sample
#  2024_Plate01   3   0     93
#  2024_Plate02   6   2     87
#  2024_Plate03   9   3     53
#  2025_Plate01  14   2     80
#  2025_Plate02   5   1     90

# =============================================================================
# next steps
#  1. if rename maps were created, run qiime feature-table rename-ids
#  2. proceed to 05_merge_plates_and_years.sh
# =============================================================================

93 + 87 + 53 + 80 + 90 = 403
3 + 6 + 9 + 14 + 5 = 37
2 + 3 + 2 + 2 + 1 = 10
#450 total
