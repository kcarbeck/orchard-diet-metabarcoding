# =============================================================================
# Identify poorly classified ASVs and prepare for BLAST
# =============================================================================
# This script:
#   1. Reads taxonomy and feature table
#   2. Identifies ASVs classified only to Arthropoda/Insecta level
#   3. Ranks them by read count
#   4. Outputs a FASTA file of top unclassified ASVs for BLAST
# =============================================================================

library(tidyverse)
library(Biostrings)  # for reading/writing FASTA

# =============================================================================
# Configuration
# =============================================================================
work_dir <- "/lustre2/home/lc736_0001/orchards/merged"
blast_dir <- file.path(work_dir, "blast_unclassified")

# How many top unclassified ASVs to BLAST?
n_to_blast <- 1000

# =============================================================================
# Load data
# =============================================================================
# Read taxonomy
tax <- read_tsv(file.path(blast_dir, "taxonomy_export/taxonomy.tsv"), 
                col_names = c("FeatureID", "Taxon", "Confidence"), 
                skip = 1, show_col_types = FALSE)

# Read feature table (skip the first line which has #OTU ID header)
ft_raw <- read_tsv(file.path(blast_dir, "feature_table.tsv"), 
                   skip = 1, show_col_types = FALSE)
names(ft_raw)[1] <- "FeatureID"

# Calculate total reads per ASV
ft <- ft_raw %>%
  mutate(total_reads = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
  select(FeatureID, total_reads)

# Read sequences
seqs <- readDNAStringSet(file.path(blast_dir, "seqs_export/dna-sequences.fasta"))

message("  Loaded ", nrow(tax), " ASVs from taxonomy")
message("  Loaded ", nrow(ft), " ASVs from feature table")
message("  Total reads in feature table: ", format(sum(ft$total_reads), big.mark = ","))

# Check how many IDs match
n_match <- sum(tax$FeatureID %in% ft$FeatureID)
message("  ", n_match, " ASVs match between taxonomy and feature table")
#  Loaded 40,294 ASVs from taxonomy
#  Loaded 33,348 ASVs from feature table
#  Total reads in feature table: 498,680,467
#  33,348 ASVs match between taxonomy and feature table


# =============================================================================
# Parse taxonomy into columns
# =============================================================================
message("\nParsing taxonomy...")

# Split taxonomy string into levels
tax_split <- tax %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", remove = FALSE) %>%
  mutate(across(Kingdom:Species, ~ str_trim(gsub("^[a-z]__", "", .x)))) %>%
  mutate(across(Kingdom:Species, ~ na_if(.x, ""))) %>%
  mutate(across(Kingdom:Species, ~ na_if(.x, "NA")))

# =============================================================================
# Identify classification depth for each ASV
# =============================================================================
message("Identifying classification depth...")

# Count how many levels are classified (not NA)
tax_split <- tax_split %>%
  rowwise() %>%
  mutate(
    depth = sum(!is.na(c(Kingdom, Phylum, Class, Order, Family, Genus, Species))),
    lowest_rank = case_when(
      !is.na(Species) ~ "Species",
      !is.na(Genus) ~ "Genus", 
      !is.na(Family) ~ "Family",
      !is.na(Order) ~ "Order",
      !is.na(Class) ~ "Class",
      !is.na(Phylum) ~ "Phylum",
      !is.na(Kingdom) ~ "Kingdom",
      TRUE ~ "Unassigned"
    )
  ) %>%
  ungroup()

# Summary of classification depth
depth_summary <- tax_split %>%
  count(lowest_rank) %>%
  arrange(match(lowest_rank, c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Unassigned")))

message("\nClassification depth summary:")
print(depth_summary)

#Classification depth summary:
# A tibble: 6 × 2
#  lowest_rank     n
#  <chr>       <int>
# 1 Species     17086
# 2 Genus        6891
# 3 Family       2557
# 4 Order        3715
# 5 Class        9238
# 6 Phylum        807

# =============================================================================
# Identify "poorly classified" ASVs
# =============================================================================
# Definition: classified only to Phylum, Class, or Order level
# (i.e., not to Family, Genus, or Species)

poorly_classified <- tax_split %>%
  filter(lowest_rank %in% c("Kingdom", "Phylum", "Class", "Order")) %>%
  left_join(ft, by = "FeatureID") %>%
  filter(!is.na(total_reads)) %>%  # keep only ASVs that are in the filtered feature table

  arrange(desc(total_reads))

message("\n", nrow(poorly_classified), " ASVs are poorly classified (Order level or higher)")
message("These account for ", format(sum(poorly_classified$total_reads), big.mark = ","), " reads (",
        round(100 * sum(poorly_classified$total_reads) / sum(ft$total_reads), 1), "% of total)")
# 10,960 ASVs are poorly classified (Order level or higher)
# These account for 127,070,149 reads (25.5% of total)

# Breakdown by lowest rank
message("\nBreakdown of poorly classified ASVs:")
poorly_classified %>%
  group_by(lowest_rank) %>%
  summarize(
    n_asvs = n(),
    total_reads = sum(total_reads),
    .groups = "drop"
  ) %>%
  arrange(desc(total_reads)) %>%
  print()
#Breakdown of poorly classified ASVs:
# A tibble: 3 × 3
#  lowest_rank n_asvs total_reads
#  <chr>        <int>       <dbl>
#1 Order         3023    72,855,347
#2 Class         7291    49,138,773
#3 Phylum         646     5,076,029

# =============================================================================
# Get top poorly classified ASVs for BLAST
# =============================================================================
top_unclassified <- poorly_classified %>%
  head(n_to_blast)

message("\nTop ", n_to_blast, " poorly classified ASVs by read count:")
print(top_unclassified %>% select(FeatureID, total_reads, lowest_rank, Order, Class, Phylum))

# =============================================================================
# Save outputs
# =============================================================================

# Save full summary
write_csv(poorly_classified, file.path(blast_dir, "poorly_classified_asvs.csv"))
message("\nSaved full list to: ", file.path(blast_dir, "poorly_classified_asvs.csv"))
#Saved full list to: /lustre2/home/lc736_0001/orchards/merged/blast_unclassified/poorly_classified_asvs.csv

# Save top ASVs for BLAST
write_csv(top_unclassified, file.path(blast_dir, "top_1000_unclassified_for_blast.csv"))

# -----------------------------------------------------------------------------
# Create FASTA file for ALL poorly classified ASVs (for command-line BLAST)
# -----------------------------------------------------------------------------
message("\nCreating FASTA for ALL poorly classified ASVs...")

all_poor_seqs <- seqs[names(seqs) %in% poorly_classified$FeatureID]

# Add read count and taxonomy info to sequence names
all_poor_names <- sapply(names(all_poor_seqs), function(id) {
  info <- poorly_classified %>% filter(FeatureID == id)
  order_info <- ifelse(is.na(info$Order), info$Class, info$Order)
  paste0(id, "|reads=", info$total_reads, "|", order_info)
})
names(all_poor_seqs) <- all_poor_names

writeXStringSet(all_poor_seqs, file.path(blast_dir, "all_poorly_classified.fasta"))
message("Saved ALL poorly classified to: ", file.path(blast_dir, "all_poorly_classified.fasta"))
message("  (", length(all_poor_seqs), " sequences)")

# -----------------------------------------------------------------------------
# Create FASTA file for top N poorly classified ASVs
# -----------------------------------------------------------------------------
top_seqs <- seqs[names(seqs) %in% top_unclassified$FeatureID]

# Add read count and taxonomy info to sequence names for easier interpretation
new_names <- sapply(names(top_seqs), function(id) {
  info <- top_unclassified %>% filter(FeatureID == id)
  order_info <- ifelse(is.na(info$Order), info$Class, info$Order)
  paste0(id, "|reads=", format(info$total_reads, big.mark = ""), "|", order_info)
})
names(top_seqs) <- new_names

writeXStringSet(top_seqs, file.path(blast_dir, "top_unclassified.fasta"))
message("Saved FASTA for BLAST to: ", file.path(blast_dir, "top_unclassified.fasta"))


# =============================================================================
# Summary statistics
# =============================================================================
message("\n============================================")
message("SUMMARY")
message("============================================")
message("")
message("Total ASVs in filtered data: ", nrow(ft))
message("Total reads: ", format(sum(ft$total_reads), big.mark = ","))
message("")
message("Poorly classified (Order or higher): ", nrow(poorly_classified), " ASVs")
message("  Reads in poorly classified: ", format(sum(poorly_classified$total_reads), big.mark = ","), 
        " (", round(100 * sum(poorly_classified$total_reads) / sum(ft$total_reads), 1), "%)")
message("")
message("Well classified (Family or lower): ", nrow(ft) - nrow(poorly_classified), " ASVs")
well_classified_reads <- sum(ft$total_reads) - sum(poorly_classified$total_reads)
message("  Reads in well classified: ", format(well_classified_reads, big.mark = ","),
        " (", round(100 * well_classified_reads / sum(ft$total_reads), 1), "%)")

# =============================================================================
# Instructions for BLAST
# =============================================================================
message("\n============================================")
message("FASTA FILES CREATED:")
message("============================================")
message("")
message("1. all_poorly_classified.fasta - ALL ", nrow(poorly_classified), " poorly classified ASVs")
message("2. top_unclassified.fasta - Top ", n_to_blast, " by read count")
message("")
message("============================================")
message("COMMAND-LINE BLAST (recommended for large files)")
message("============================================")
message("")
message("# For ALL poorly classified sequences:")
message("cd ", blast_dir)
message("")
message("# Using remote BLAST (no local DB needed, but slower):")
message("blastn -query all_poorly_classified.fasta \\")
message("       -db nt -remote \\")
message("       -max_target_seqs 3 \\")
message("       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames' \\")
message("       -out blast_all_results.tsv")
message("")
message("# Or if you have a local nt database:")
message("blastn -query all_poorly_classified.fasta \\")
message("       -db /path/to/nt \\")
message("       -num_threads 8 \\")
message("       -max_target_seqs 3 \\")
message("       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames' \\")
message("       -out blast_all_results.tsv")
message("")
message("============================================")
message("")
message("Common reasons for poor classification:")
message("  - Species not in reference database (regional/rare species)")
message("  - Sequence divergent from references")
message("  - Under-represented groups (some mite families, springtails)")
message("  - Chimeric sequences (should be rare after DADA2)")
message("")




# ============================================
# SUMMARY
# ============================================
# 
# Total ASVs in filtered data: 33348
# Total reads: 498,680,467
# 
# Poorly classified (Order or higher): 10960 ASVs
#   Reads in poorly classified: 127,070,149 (25.5%)
# 
# Well classified (Family or lower): 22388 ASVs
#   Reads in well classified: 371,610,318 (74.5%)