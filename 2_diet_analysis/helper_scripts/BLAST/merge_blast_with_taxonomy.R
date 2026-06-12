# =============================================================================
# Merge BLAST results with taxonomy
# =============================================================================
# This script:
#   1. Reads BLAST results for poorly classified ASVs
#   2. Filters for high-quality hits (pident, coverage)
#   3. Looks up full taxonomy from NCBI using taxids (parallelized)
#   4. Merges with original classifier taxonomy
#   5. Creates enhanced taxonomy for downstream analysis
# =============================================================================
#module load R/4.2.3
library(tidyverse)
library(taxize)      # for NCBI taxonomy lookups
library(furrr)       # for parallel processing
library(future)      # for parallel backend

# =============================================================================
# Configuration
# =============================================================================
work_dir <- "/lustre2/home/lc736_0001/orchards/merged"
blast_dir <- file.path(work_dir, "blast_unclassified")

# Quality thresholds for accepting BLAST hits
MIN_LENGTH <- 150     # minimum alignment length (bp) - ~83% of amplicon
MIN_COVERAGE <- 0.83  # minimum query coverage (150/180bp)
QUERY_LENGTH <- 180   # actual amplicon length for this study

# -----------------------------------------------------------------------------
# Percent identity thresholds for taxonomic rank assignment
# -----------------------------------------------------------------------------
# Based on standard cutoffs for COI/metabarcoding markers:
#   Species: 97-99%+ identity
#   Genus:   94-97%  identity  
#   Family:  90-94%  identity
#
# We use the LOWER bound of each range to assign that rank:
#   >= 97% -> species level
#   >= 94% -> genus level (but not species)
#   >= 90% -> family level (but not genus/species)
#   <  90% -> reject (not reliable enough)
# -----------------------------------------------------------------------------
PIDENT_SPECIES <- 97   # minimum for species-level assignment
PIDENT_GENUS   <- 94   # minimum for genus-level assignment
PIDENT_FAMILY  <- 90   # minimum for family-level assignment (also minimum to keep hit)

# Parallelization settings for NCBI lookups
N_CORES <- 20          # number of cores to use (adjust for your cluster)
BATCH_SIZE <- 50      # number of taxids to query at once (NCBI limit-friendly)

# =============================================================================
# Load data
# =============================================================================
message("Loading data...")

# BLAST results
# Column names based on outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames'
blast_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                "qstart", "qend", "sstart", "send", "evalue", "bitscore", 
                "staxids", "sscinames")

blast <- read_tsv(file.path(blast_dir, "blast_all_results.tsv"), 
                  col_names = blast_cols, show_col_types = FALSE)

message("  Loaded ", nrow(blast), " BLAST hits for ", n_distinct(blast$qseqid), " ASVs")
#  Loaded 49777 BLAST hits for 9699 ASVs

# Original taxonomy
tax_orig <- read_tsv(file.path(blast_dir, "taxonomy_export/taxonomy.tsv"),
                     col_names = c("FeatureID", "Taxon", "Confidence"),
                     skip = 1, show_col_types = FALSE)

# Poorly classified ASVs (from previous script)
poorly_classified <- read_csv(file.path(blast_dir, "poorly_classified_asvs.csv"),
                              show_col_types = FALSE)

message("  ", nrow(poorly_classified), " poorly classified ASVs")
#  10960 poorly classified ASVs

# =============================================================================
# Parse BLAST query IDs
# =============================================================================
# Query IDs have format: FeatureID|reads=XXXX|Order
# Extract the actual FeatureID

blast <- blast %>%
  mutate(
    FeatureID = str_extract(qseqid, "^[^|]+"),
    query_reads = as.numeric(str_extract(qseqid, "(?<=reads=)[0-9]+"))
  )

# =============================================================================
# Filter BLAST hits
# =============================================================================
message("\nFiltering BLAST hits...")
message("  Minimum thresholds: pident >= ", PIDENT_FAMILY, "%, length >= ", MIN_LENGTH, 
        "bp, coverage >= ", MIN_COVERAGE * 100, "%")
message("  Tiered taxonomy assignment:")
message("    >= ", PIDENT_SPECIES, "% -> species level")
message("    >= ", PIDENT_GENUS, "% -> genus level")
message("    >= ", PIDENT_FAMILY, "% -> family level")

blast_filtered <- blast %>%
  mutate(coverage = length / QUERY_LENGTH) %>%
  filter(
    pident >= PIDENT_FAMILY,  # use family threshold as minimum to keep
    length >= MIN_LENGTH,
    coverage >= MIN_COVERAGE
  )

message("  ", nrow(blast_filtered), " hits passed filters")
message("  ", n_distinct(blast_filtered$FeatureID), " ASVs have good BLAST hits")

# Get best hit per ASV (highest bitscore)
blast_best <- blast_filtered %>%
  group_by(FeatureID) %>%
  slice_max(bitscore, n = 1, with_ties = FALSE) %>%
  ungroup()

message("  ", nrow(blast_best), " ASVs with best hits selected")
#  23499 hits passed filters
#  5149 ASVs have good BLAST hits
#  5149 ASVs with best hits selected

# =============================================================================
# Look up full taxonomy from NCBI using taxids (PARALLELIZED)
# =============================================================================
message("\n============================================")
message("NCBI TAXONOMY LOOKUP")
message("============================================")

# Check if we have a cached version (to avoid re-querying NCBI)
ncbi_cache_file <- file.path(blast_dir, "ncbi_taxonomy_cache.csv")

# Get unique taxids from best hits
# Handle multiple taxids per hit (some have "taxid1;taxid2" format)
unique_taxids <- blast_best %>%
  mutate(taxid_clean = str_extract(staxids, "^[0-9]+")) %>%  # take first taxid if multiple
  pull(taxid_clean) %>%
  unique() %>%
  na.omit() %>%
  as.integer()

message("\nFound ", length(unique_taxids), " unique taxids to look up")
# Found 556 unique taxids to look up

# Check if cache exists and has all our taxids
use_cache <- FALSE
if (file.exists(ncbi_cache_file)) {
  ncbi_cache <- read_csv(ncbi_cache_file, show_col_types = FALSE)
  cached_taxids <- ncbi_cache$taxid
  missing_taxids <- setdiff(unique_taxids, cached_taxids)
  
  if (length(missing_taxids) == 0) {
    message("  Using cached NCBI taxonomy (", nrow(ncbi_cache), " taxids)")
    message("  To re-query NCBI, delete: ", ncbi_cache_file)
    ncbi_taxonomy <- ncbi_cache
    use_cache <- TRUE
  } else {
    message("  Cache exists but missing ", length(missing_taxids), " taxids - will query all")
  }
}

if (!use_cache) {
  message("Using ", N_CORES, " cores for parallel processing")
  
  # Set up parallel backend
  plan(multisession, workers = N_CORES)
  
  # Function to safely get classification for a batch of taxids
  get_taxonomy_batch <- function(taxid_batch) {
    tryCatch({
      # Query NCBI for classification
      result <- classification(taxid_batch, db = "ncbi")
      
      # Parse results into a tidy dataframe
      map_dfr(names(result), function(tid) {
        cl <- result[[tid]]
        if (is.null(cl) || inherits(cl, "logical") || nrow(cl) == 0) {
          return(tibble(taxid = as.integer(tid)))  # return empty row if failed
        }
        
        # Extract each rank we care about
        tibble(
          taxid = as.integer(tid),
          ncbi_kingdom = cl$name[cl$rank == "kingdom"][1],
          ncbi_phylum = cl$name[cl$rank == "phylum"][1],
          ncbi_class = cl$name[cl$rank == "class"][1],
          ncbi_order = cl$name[cl$rank == "order"][1],
          ncbi_family = cl$name[cl$rank == "family"][1],
          ncbi_genus = cl$name[cl$rank == "genus"][1],
          ncbi_species = cl$name[cl$rank == "species"][1]
        )
      })
    }, error = function(e) {
      message("  Error in batch: ", e$message)
      tibble(taxid = taxid_batch)  # return taxids with no taxonomy on error
    })
  }
  
  # Split taxids into batches
  taxid_batches <- split(unique_taxids, ceiling(seq_along(unique_taxids) / BATCH_SIZE))
  message("Split into ", length(taxid_batches), " batches of ~", BATCH_SIZE, " taxids each")
  
  # Process batches in parallel with progress
  message("\nQuerying NCBI (this may take a few minutes)...")
  start_time <- Sys.time()
  
  ncbi_taxonomy <- future_map_dfr(
    taxid_batches, 
    get_taxonomy_batch,
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  )
  
  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
  message("  Completed in ", elapsed, " minutes")
  message("  Retrieved taxonomy for ", sum(!is.na(ncbi_taxonomy$ncbi_family)), " / ", 
          length(unique_taxids), " taxids")
  
  # Reset to sequential processing
  plan(sequential)
  
  # Save the NCBI lookup results for future use (avoids re-querying)
  write_csv(ncbi_taxonomy, ncbi_cache_file)
  message("  Cached NCBI results to: ", ncbi_cache_file)
}  # end if (!use_cache)
#   Retrieved taxonomy for 532 / 556 taxids

# Join NCBI taxonomy to blast_best
blast_best <- blast_best %>%
  mutate(taxid_clean = as.integer(str_extract(staxids, "^[0-9]+"))) %>%
  left_join(ncbi_taxonomy, by = c("taxid_clean" = "taxid"))

# =============================================================================
# Parse species names from BLAST
# =============================================================================
# sscinames might be empty if taxdb wasn't available
# We can try to extract genus/species from sseqid if needed

blast_best <- blast_best %>%
  mutate(
    # Try to clean up species name
    blast_species = case_when(
      !is.na(sscinames) & sscinames != "" ~ sscinames,
      TRUE ~ NA_character_
    ),
    # Extract genus (first word of species name)
    blast_genus = str_extract(blast_species, "^[A-Z][a-z]+"),
    # Clean up - remove subspecies/strain info
    blast_species_clean = str_extract(blast_species, "^[A-Z][a-z]+ [a-z]+"),
    # Use NCBI genus/species if available (more reliable)
    blast_genus = coalesce(ncbi_genus, blast_genus),
    blast_species_clean = coalesce(ncbi_species, blast_species_clean)
  )

# =============================================================================
# Summary of BLAST results
# =============================================================================
message("\n============================================")
message("BLAST RESULTS SUMMARY")
message("============================================")

# How many ASVs got identified?
n_identified <- sum(!is.na(blast_best$blast_species))
message("")
message("ASVs with species-level BLAST ID: ", n_identified, " / ", nrow(poorly_classified))
message("ASVs without good BLAST hit: ", nrow(poorly_classified) - nrow(blast_best))

# Top species identified
message("\nTop 20 species identified by BLAST:")
blast_best %>%
  filter(!is.na(blast_species_clean)) %>%
  count(blast_species_clean, sort = TRUE) %>%
  head(20) %>%
  print()

#  1 Rhagio tringarius            345
#  2 Barypeithes pellucidus       240
#  3 Trachelipus rathkii          204
#  4 Sylvicola alternatus         198
#  5 Amphipoea velata             196
#  6 Tarsonemidae sp              175
#  7 Sphaerosorus coelastroides   149
#  8 Cecidomyiidae sp             144
#  9 Nannochloris sp              142
# 10 Philoscia muscorum           131
# 11 Cylindroiulus punctatus       92
# 12 Porcellio spinicornis         85
# 13 Oribatella sp                 76
# 14 Valenzuela flavidus           67
# 15 Steneotarsonemus laticeps     65
# 16 Teliapsocus sp                64
# 17 Sialis sp                     63
# 18 Adineta vaga                  62
# 19 Chrysopilus asiliformis       59
# 20 Peloptulus phaenotus          51

# =============================================================================
# Merge with original taxonomy
# =============================================================================
message("\nMerging BLAST results with original taxonomy...")

# Parse original taxonomy into columns
tax_parsed <- tax_orig %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", remove = FALSE) %>%
  mutate(across(Kingdom:Species, ~ str_trim(gsub("^[a-z]__", "", .x)))) %>%
  mutate(across(Kingdom:Species, ~ na_if(.x, ""))) %>%
  mutate(across(Kingdom:Species, ~ na_if(.x, "NA")))

# Add BLAST results with TIERED taxonomy assignment based on percent identity
# -----------------------------------------------------------------------------
# Logic:
#   - When BLAST has a hit, use NCBI taxonomy COMPLETELY for higher ranks
#     (Kingdom through Family) to avoid chimeric taxonomy
#   - For Genus: only use BLAST if pident >= PIDENT_GENUS (94%)
#   - For Species: only use BLAST if pident >= PIDENT_SPECIES (97%)
#   - If no BLAST hit, keep original classifier taxonomy
#
# This prevents over-confident species assignments from low-identity matches
# while still using BLAST to correct misidentified organisms at higher ranks
# -----------------------------------------------------------------------------
tax_enhanced <- tax_parsed %>%
  left_join(
    blast_best %>% select(FeatureID,
                          blast_genus, blast_species_clean,
                          ncbi_kingdom, ncbi_phylum, ncbi_class, ncbi_order, ncbi_family,
                          pident, bitscore),
    by = "FeatureID"
  ) %>%
  mutate(
    # Did BLAST provide a hit? (blast_best is already filtered to >= PIDENT_FAMILY)
    blast_has_hit = !is.na(pident),
    # Determine what taxonomic level is supported by the pident
    blast_supports_species = blast_has_hit & pident >= PIDENT_SPECIES,
    blast_supports_genus = blast_has_hit & pident >= PIDENT_GENUS,
    # KEY CHECK: Does BLAST confirm the same organism as the classifier?
    # Compare at family level (or phylum if family not available)
    # If they match, we can trust the classifier's more detailed assignments
    # If they differ, the classifier was wrong and we should not keep genus/species
    blast_confirms_classifier = case_when(
      !blast_has_hit ~ FALSE,
      !is.na(ncbi_family) & !is.na(Family) ~ ncbi_family == Family,
      !is.na(ncbi_phylum) & !is.na(Phylum) ~ ncbi_phylum == Phylum,
      TRUE ~ FALSE
    ),
    # Higher ranks (Kingdom through Family): use NCBI when BLAST has any hit
    # This ensures we don't mix "Fungi kingdom" with "Arthropoda phylum"
    Kingdom_enhanced = case_when(
      blast_has_hit & !is.na(ncbi_kingdom) ~ ncbi_kingdom,
      TRUE ~ Kingdom
    ),
    Phylum_enhanced = case_when(
      blast_has_hit & !is.na(ncbi_phylum) ~ ncbi_phylum,
      TRUE ~ Phylum
    ),
    Class_enhanced = case_when(
      blast_has_hit & !is.na(ncbi_class) ~ ncbi_class,
      TRUE ~ Class
    ),
    Order_enhanced = case_when(
      blast_has_hit & !is.na(ncbi_order) ~ ncbi_order,
      TRUE ~ Order
    ),
    Family_enhanced = case_when(
      blast_has_hit & !is.na(ncbi_family) ~ ncbi_family,
      TRUE ~ Family
    ),
    # GENUS/SPECIES LOGIC:
    # 1. If BLAST pident >= threshold: use BLAST assignment
    # 2. If BLAST confirms classifier (same family): keep classifier assignment
    # 3. If BLAST contradicts classifier: blank it (classifier was wrong)
    # 4. If no BLAST hit: keep classifier assignment
    # Genus: use BLAST if >= 94%, else keep classifier if BLAST confirms same family
    Genus_enhanced = case_when(
      blast_supports_genus & !is.na(blast_genus) ~ blast_genus,
      blast_confirms_classifier & !is.na(Genus) ~ Genus,
      blast_has_hit ~ NA_character_,
      TRUE ~ Genus
    ),

    # Species: use BLAST if >= 97%, else keep classifier if BLAST confirms same family
    Species_enhanced = case_when(
      blast_supports_species & !is.na(blast_species_clean) ~ blast_species_clean,
      blast_confirms_classifier & !is.na(Species) ~ Species,
      blast_has_hit ~ NA_character_,
      TRUE ~ Species
    ),

    # Track source and depth of classification
    classification_source = case_when(
      blast_supports_species ~ "BLAST_species",
      blast_supports_genus ~ "BLAST_genus",
      blast_has_hit & blast_confirms_classifier & !is.na(Species) ~ "classifier_confirmed",
      blast_has_hit & blast_confirms_classifier & !is.na(Genus) ~ "classifier_confirmed",
      blast_has_hit ~ "BLAST_family",
      !is.na(Species) ~ "classifier",
      !is.na(Genus) ~ "classifier",
      TRUE ~ "unclassified"
    ),
    # Store the pident for reference
    blast_pident = pident
  ) %>%
  select(-blast_has_hit, -blast_supports_species, -blast_supports_genus, -blast_confirms_classifier)

# Summary of enhancement
message("\n============================================")
message("TIERED TAXONOMY ASSIGNMENT SUMMARY")
message("============================================")
message("")
message("Classification source breakdown:")
tax_enhanced %>%
  count(classification_source) %>%
  arrange(desc(n)) %>%
  print()

# Breakdown of BLAST assignments by pident
message("")
message("BLAST assignments by percent identity:")
blast_summary <- tax_enhanced %>%
  filter(str_detect(classification_source, "BLAST|classifier_confirmed")) %>%
  summarize(
    total_blast_hits = n(),
    species_level = sum(classification_source == "BLAST_species"),
    genus_level = sum(classification_source == "BLAST_genus"),
    classifier_confirmed = sum(classification_source == "classifier_confirmed"),
    family_level = sum(classification_source == "BLAST_family"),
    mean_pident = round(mean(blast_pident, na.rm = TRUE), 1),
    min_pident = round(min(blast_pident, na.rm = TRUE), 1),
    max_pident = round(max(blast_pident, na.rm = TRUE), 1)
  )
message("  Total ASVs with BLAST hits: ", blast_summary$total_blast_hits)
message("    BLAST species-level (>=", PIDENT_SPECIES, "%): ", blast_summary$species_level)
message("    BLAST genus-level (>=", PIDENT_GENUS, "%):     ", blast_summary$genus_level)
message("    Classifier confirmed*:        ", blast_summary$classifier_confirmed)
message("    BLAST family-level (>=", PIDENT_FAMILY, "%):   ", blast_summary$family_level)
message("  Percent identity range: ", blast_summary$min_pident, "% - ", blast_summary$max_pident, "%")
message("  Mean percent identity: ", blast_summary$mean_pident, "%")
message("")
message("  *Classifier confirmed: BLAST hit <", PIDENT_GENUS, "% but confirmed same family,")
message("   so kept classifier's more detailed genus/species assignment")

# Classification source breakdown:
# # A tibble: 5 × 2
#   classification_source     n
#   <chr>                 <int>
# 1 classifier            23977
# 2 unclassified          11168
# 3 BLAST_species          3726
# 4 BLAST_genus             734
# 5 BLAST_family            689
# 
# BLAST assignments by percent identity:
#   Total ASVs with BLAST hits: 5149
#     BLAST species-level (>=97%): 3726
#     BLAST genus-level (>=94%):     734
#     Classifier confirmed*:        0
#     BLAST family-level (>=90%):   689
#   Percent identity range: 90% - 100%
#   Mean percent identity: 97.7%
# 
#   *Classifier confirmed: BLAST hit <94% but confirmed same family, so kept classifier's more detailed genus/species assignment


# -----------------------------------------------------------------------------
# QC: Check for cases where BLAST identified a DIFFERENT organism than classifier
# -----------------------------------------------------------------------------
# This happens when the classifier misidentified the sequence (e.g., fungus as Arthropoda)
# With the fixed logic, we now use BLAST taxonomy fully in these cases

# Join original classifier phylum with BLAST phylum to compare
qc_mismatch <- tax_enhanced %>%
  filter(str_detect(classification_source, "BLAST")) %>%
  filter(!is.na(Phylum) & !is.na(ncbi_phylum)) %>%
  filter(Phylum != ncbi_phylum)

if (nrow(qc_mismatch) > 0) {
  message("\n============================================")
  message("QC: BLAST corrected ", nrow(qc_mismatch), " ASVs with WRONG phylum from classifier")
  message("============================================")
  message("These were misidentified by the classifier (e.g., fungi classified as Arthropoda)")
  message("BLAST identified them correctly - using BLAST taxonomy fully for these ASVs")
  message("")
  
  # Show breakdown of corrections
  mismatch_summary <- qc_mismatch %>%
    count(Phylum, ncbi_phylum, name = "n_asvs") %>%
    arrange(desc(n_asvs))
  
  message("Corrections made (classifier phylum -> BLAST phylum):")
  print(mismatch_summary %>% head(20))
}

# Corrections made (classifier phylum -> BLAST phylum):
# # A tibble: 17 × 3
#    Phylum     ncbi_phylum     n_asvs
#    <chr>      <chr>            <int>
#  1 Arthropoda Chlorophyta        226
#  2 Arthropoda Rotifera           162
#  3 Arthropoda Mollusca           124
#  4 Arthropoda Cercozoa            81
#  5 Arthropoda Discosea            66
#  6 Arthropoda Annelida            55
#  7 Arthropoda Chordata            50
#  8 Arthropoda Oomycota            30
#  9 Arthropoda Streptophyta        16
# 10 Arthropoda Acanthocephala       6
# 11 Arthropoda Ascomycota           6
# 12 Arthropoda Nematoda             6
# 13 Arthropoda Pseudomonadota       5
# 14 Arthropoda Bacillariophyta      1
# 15 Arthropoda Bacillota            1
# 16 Arthropoda Chlamydiota          1
# 17 Arthropoda Nemertea             1

# =============================================================================
# Create enhanced taxonomy string
# =============================================================================
# Format: k__Kingdom;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species
# Uses ALL enhanced levels (filled in from BLAST where classifier was missing)
# Empty string for missing ranks (not "NA") to avoid downstream parsing issues

tax_enhanced <- tax_enhanced %>%
  mutate(
    Taxon_enhanced = paste(
      paste0("k__", ifelse(is.na(Kingdom_enhanced), "", Kingdom_enhanced)),
      paste0("p__", ifelse(is.na(Phylum_enhanced), "", Phylum_enhanced)),
      paste0("c__", ifelse(is.na(Class_enhanced), "", Class_enhanced)),
      paste0("o__", ifelse(is.na(Order_enhanced), "", Order_enhanced)),
      paste0("f__", ifelse(is.na(Family_enhanced), "", Family_enhanced)),
      paste0("g__", ifelse(is.na(Genus_enhanced), "", Genus_enhanced)),
      paste0("s__", ifelse(is.na(Species_enhanced), "", Species_enhanced)),
      sep = ";"
    )
  )

# =============================================================================
# Save outputs
# =============================================================================
output_dir <- file.path(blast_dir, "enhanced_taxonomy")
dir.create(output_dir, showWarnings = FALSE)

# Full enhanced taxonomy table
write_csv(tax_enhanced, file.path(output_dir, "taxonomy_enhanced_full.csv"))

# QIIME2-compatible taxonomy file (just Feature ID, Taxon, Confidence)
# Note: QIIME2 requires "Feature ID" with a space, not "FeatureID"
tax_for_qiime <- tax_enhanced %>%
  select(`Feature ID` = FeatureID, Taxon = Taxon_enhanced, Confidence) %>%
  mutate(Confidence = ifelse(is.na(Confidence), 0.5, Confidence))  # default confidence for BLAST

write_tsv(tax_for_qiime, file.path(output_dir, "taxonomy_enhanced.tsv"))

# Summary of BLAST hits
write_csv(blast_best, file.path(output_dir, "blast_best_hits.csv"))

message("\n============================================")
message("FILES SAVED:")
message("============================================")
message("")
message("1. taxonomy_enhanced_full.csv - Full table with all columns")
message("2. taxonomy_enhanced.tsv - QIIME2-compatible format")
message("3. blast_best_hits.csv - Best BLAST hit per ASV")
message("")
message("Location: ", output_dir)

#============================================
#FILES SAVED:
#============================================
#1. taxonomy_enhanced_full.csv - Full table with all columns
#2. taxonomy_enhanced.tsv - QIIME2-compatible format
#3. blast_best_hits.csv - Best BLAST hit per ASV
#Location: /lustre2/home/lc736_0001/orchards/merged/blast_unclassified/enhanced_taxonomy

# =============================================================================
# Summary statistics
# =============================================================================
message("\n============================================")
message("FINAL SUMMARY")
message("============================================")
message("")
message("Tiered thresholds used:")
message("  Species: >= ", PIDENT_SPECIES, "% identity")
message("  Genus:   >= ", PIDENT_GENUS, "% identity")
message("  Family:  >= ", PIDENT_FAMILY, "% identity")

# Count classifications at each level
count_at_level <- function(df, col) {
  sum(!is.na(df[[col]]) & df[[col]] != "NA")
}

message("")
message("Original classifier taxonomy:")
message("  Kingdom-level: ", count_at_level(tax_parsed, "Kingdom"))
message("  Phylum-level:  ", count_at_level(tax_parsed, "Phylum"))
message("  Class-level:   ", count_at_level(tax_parsed, "Class"))
message("  Order-level:   ", count_at_level(tax_parsed, "Order"))
message("  Family-level:  ", count_at_level(tax_parsed, "Family"))
message("  Genus-level:   ", count_at_level(tax_parsed, "Genus"))
message("  Species-level: ", count_at_level(tax_parsed, "Species"))

message("")
message("After adding BLAST + NCBI taxonomy (with tiered thresholds):")
message("  Kingdom-level: ", count_at_level(tax_enhanced, "Kingdom_enhanced"))
message("  Phylum-level:  ", count_at_level(tax_enhanced, "Phylum_enhanced"))
message("  Class-level:   ", count_at_level(tax_enhanced, "Class_enhanced"))
message("  Order-level:   ", count_at_level(tax_enhanced, "Order_enhanced"))
message("  Family-level:  ", count_at_level(tax_enhanced, "Family_enhanced"))
message("  Genus-level:   ", count_at_level(tax_enhanced, "Genus_enhanced"), 
        " (BLAST only if >= ", PIDENT_GENUS, "%)")
message("  Species-level: ", count_at_level(tax_enhanced, "Species_enhanced"),
        " (BLAST only if >= ", PIDENT_SPECIES, "%)")

message("")
message("Changes from BLAST + NCBI (may be negative at genus/species due to tiered thresholds):")
message("  ", sprintf("%+d", count_at_level(tax_enhanced, "Kingdom_enhanced") - count_at_level(tax_parsed, "Kingdom")), " ASVs with Kingdom")
message("  ", sprintf("%+d", count_at_level(tax_enhanced, "Phylum_enhanced") - count_at_level(tax_parsed, "Phylum")), " ASVs with Phylum")
message("  ", sprintf("%+d", count_at_level(tax_enhanced, "Class_enhanced") - count_at_level(tax_parsed, "Class")), " ASVs with Class")
message("  ", sprintf("%+d", count_at_level(tax_enhanced, "Order_enhanced") - count_at_level(tax_parsed, "Order")), " ASVs with Order")
message("  ", sprintf("%+d", count_at_level(tax_enhanced, "Family_enhanced") - count_at_level(tax_parsed, "Family")), " ASVs with Family")
message("  ", sprintf("%+d", count_at_level(tax_enhanced, "Genus_enhanced") - count_at_level(tax_parsed, "Genus")), " ASVs with Genus")
message("  ", sprintf("%+d", count_at_level(tax_enhanced, "Species_enhanced") - count_at_level(tax_parsed, "Species")), " ASVs with Species")

# ============================================
# FINAL SUMMARY
# ============================================
# 
# Tiered thresholds used:
#   Species: >= 97% identity
#   Genus:   >= 94% identity
#   Family:  >= 90% identity
# 
# Original classifier taxonomy:
#   Kingdom-level: 0
#   Phylum-level:  40294
#   Class-level:   39487
#   Order-level:   30249
#   Family-level:  26534
#   Genus-level:   23977
#   Species-level: 17086
# 
# After adding BLAST + NCBI taxonomy (with tiered thresholds):
#   Kingdom-level: 4752
#   Phylum-level:  40294
#   Class-level:   39724
#   Order-level:   32896
#   Family-level:  31617
#   Genus-level:   28437 (BLAST only if >= 94%)
#   Species-level: 20812 (BLAST only if >= 97%)
# 
# Changes from BLAST + NCBI (may be negative at genus/species due to tiered thresholds):
#   +4752 ASVs with Kingdom
#   +0 ASVs with Phylum
#   +237 ASVs with Class
#   +2647 ASVs with Order
#   +5083 ASVs with Family
#   +4460 ASVs with Genus
#   +3726 ASVs with Species






# =============================================================================
# Next steps
# =============================================================================
message("\n============================================")
message("NEXT STEPS:")
message("============================================")
message("")
message("To use enhanced taxonomy in QIIME2:")
message("  1. Import the enhanced taxonomy:")
message("     qiime tools import \\")
message("       --type 'FeatureData[Taxonomy]' \\")
message("       --input-format TSVTaxonomyFormat \\")
message("       --input-path ", file.path(output_dir, "taxonomy_enhanced.tsv"), " \\")
message("       --output-path taxonomy_enhanced.qza")
message("")
message("  2. Use taxonomy_enhanced.qza instead of taxonomy_merged.qza")
message("     in your downstream analyses")
message("")

