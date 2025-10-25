# subset databases for NE states using GBIF
# katherine carbeck
# 2 oct 2025

# this uses the gbif_filter_states.R script in the helper_scripts directory to filter for only records that occur within the range of your study area. It retains database sequences only if the taxon occurs in GBIF within the selected US state(s) and year range (has coordinates). For each record, matching is attempted at the most specific available rank, in order: species → genus → family → order → class → phylum. A record falls back to a broader rank only when more-specific ranks are missing (blank/NA). If a more-specific rank is present but not found in GBIF, the record is excluded and does not fall back. Species-level matches require a binomial name. “Presence in GBIF” reflects any occurrence record in GBIF under the filters; it does not imply nativeness or abundance. Using GBIF bulk download ensures full coverage of matching occurrences; API mode may paginate and is not recommended when completeness is critical.
GBIF_PWD_FILE="/lustre2/home/lc736_0001/diet/songbird_coi_database/.secrets/gbif_pwd.txt"
[ -f "$GBIF_PWD_FILE" ] || { echo "Missing $GBIF_PWD_FILE" >&2; exit 1; }
export GBIF_PWD="$(tr -d '\r\n' < "$GBIF_PWD_FILE")"

##############################################################################
#*                   TEST SUBSET OF DATABASE
##############################################################################
# Test with first 1000 lines
head -1000 orchard_database_011025.txt > test_input.txt

# Run a quick diagnostic to see how many records GBIF has
Rscript -e "
library(rgbif)
count <- occ_count(country='US', stateProvince='New York', year='2020,2025', hasCoordinate=TRUE)
cat('GBIF has', format(count, big.mark=','), 'occurrences for NY 2020-2025\n')
"
#GBIF has 36,699,530 occurrences for NY 2020-2025

# Run optimized script  
time Rscript gbif_filter_states_new.R \
  --input test_input.txt \
  --output test_output.txt \
  --states "NY" \
  --year-start 2020 \
  --year-end 2025 \
  --threads 12 \
  --download always \
  --gbif-user "jenwalsh123" \
  --gbif-pwd "$GBIF_PWD" \
  --gbif-email "jlw395@cornell.edu"


##############################################################################
#*              GBIF FILTER FOR FULL DATABASE
##############################################################################
# Use 'time' and 'nohup' for long runs
nohup time Rscript gbif_filter_states.R \
  --input orchard_database_011025.txt \
  --output orchard_database_NE_states_021025.txt \
  --states "NY,PA,MA,VT,NH" \
  --year-start 1970 \
  --year-end 2025 \
  --threads 12 \
  --cache-dir .gbif_cache \
  --gbif-user "jenwalsh123" \
  --gbif-pwd "$GBIF_PWD" \
  --gbif-email "jlw395@cornell.edu" \
  --download always \
  > gbif_filter.log 2>&1 &
#[2025-10-02 15:14:54] Target states: New York, Pennsylvania, Massachusetts, Vermont, New Hampshire
#[2025-10-02 15:14:54] Using 12 threads for parallel processing
#[2025-10-02 15:14:54] Reading database: orchard_database_011025.txt
#[2025-10-02 15:15:07] Loaded 2677174 records with 11 columns
#[2025-10-02 15:15:11] Total records: 2677174; with binomial species: 966735
#[2025-10-02 15:15:14] Unique species in database: 431327
#[2025-10-02 15:15:14] Requesting GBIF bulk download via occ_download()...
#[2025-10-02 15:15:18] GBIF download key: 0041022-250920141307145
#status: preparing
#status: running
#status: succeeded
# status: succeeded
# download is done, status: succeeded
# Download file size: 21147.1 MB
# On disk at ./0041022-250920141307145.zip
# [2025-10-02 15:49:30] Importing taxonomy columns from GBIF download...
# [2025-10-02 16:07:13] GBIF unique species across states (download): 33517
# Warning message:
# In data.table::fread(cmd = paste("unzip -p", shQuote(zip_path),  :
#   Found and resolved improper quoting out-of-sample. First healed line 1246: <<2573123739 e42c4be9-dc6d-466b-8df2-60ac8c47fadd     # eb9451db-c8c1-49c7-9769-0d5b03d16f82    Plantae Tracheophyta       Magnoliopsida   Asterales       Campanulaceae   Campanula       # Campanula rotundifolia             SPECIES Campanula rotundifolia L.       Campanula rotundifolia  L.US       "The cut", 1.5 mi SE along # R.R. tracks from East Allen St.      Vermont PRESENT   fee1eb85-53bb-4edb-91cc-de55f11b3cee     44.49492        -73.18447       1538.# 0            1991-09-21       21      9       1991    5410907 5410907 PRESERVED_SPECIMEN      VT        UVMVT097626      # 47                      CC_BY_4_0       University of Vermont   M. H. >>. If the fields are not quoted (e.g. field separator does not # appear within any field), try quote="" to avoid this warning.
# [2025-10-02 16:07:13] Species in DB with GBIF occurrences in target states: 8239
# [2025-10-02 16:07:15] Filtered records: 1605261; by-species: 93299; by-genus: 234864; by-family: 1179583; by-order: 97512; by-class: 3; # by-phylum: 0
# [2025-10-02 16:07:15] Wrote filtered database to: orchard_database_NE_states_021025.txt
# [2025-10-02 16:07:15] Complete: 1605261 of 2677174 records retained
# 917.23user 118.02system 52:23.70elapsed 32%CPU (0avgtext+0avgdata 100545108maxresident)k
# 145909728inputs+224112456outputs (26982major+7743077minor)pagefaults 0swaps


# monitor progress
tail -f gbif_filter.log

# poll status (safe to run anytime)
Rscript -e 'library(rgbif); m<-occ_download_meta("0041022-250920141307145"); cat(m$status, "\n")'
#   RUNNING

# check cache size
du -sh .gbif_cache 2>/dev/null
#    24K	.gbif_cache

#if status switches to SUCCEEDED but the log doesn’t advance:
# look for the downloaded zip that rgbif fetched
ls -lh *.zip 2>/dev/null || true

# or show more detail on the download record
Rscript -e 'library(rgbif); m<-occ_download_meta("0041022-250920141307145"); str(m)'

##############################################################################
#*                   CHECK OUTPUT
##############################################################################
# confirm row counts (original vs filtered)
wc -l orchard_database_011025.txt orchard_database_NE_states_021025.txt
# 2677174 orchard_database_011025.txt
# 1605262 orchard_database_NE_states_021025.txt

# # how many unique species GBIF returned (should match the log)
awk -F'\t' '{print $10}' orchard_database_NE_states_021025.txt \
  | grep -E '^[^[:space:]]+ [^[:space:]]+$' | sort -u | wc -l
# 8239

# see breakdown of how many rows came from each taxonomic level 
grep 'Filtered records:' -n gbif_filter.log
#22:[2025-10-02 16:07:15] Filtered records: 1605261; by-species: 93299; by-genus: 234864; by-family: 1179583; by-order: 97512; by-class: 3; by-phylum: 0

# record provenance
printf "GBIF key\t0041022-250920141307145\nstates\tNY,PA,MA,VT,NH\nyears\t1970–2025\n" > orchard_database_NE_states_021025.META


##############################################################################
#*                  REFORMAT OUTPUT -__-
##############################################################################
# the output doesn't match the original format as expected, so now we need to reformat slightly to match expectations
# reads the headered/quoted output, converts all blank strings to NA,writes a TSV with no header, no quotes, and literal NA, preserves the exact 11-column order

Rscript -e 'library(data.table);
x <- fread("orchard_database_NE_states_021025.txt", sep="\t", header=TRUE, quote="\"");
# convert all "" to NA so fwrite can serialize them as "NA"
for (j in seq_along(x)) {
  if (is.character(x[[j]])) set(x, which(x[[j]] == ""), j, NA_character_)
}
# enforce the expected column order (11 cols)
setcolorder(x, c("accession","provided_name","taxid","note","phylum","class","order","family","genus","species","sequence"))
# write EXACT CRABS-style: tsv, no header, no quotes, NA as literal
fwrite(x, file="orchard_database_NE_states_021025_formatted.txt",
       sep="\t", quote=FALSE, na="NA", col.names=FALSE, showProgress=FALSE)'


# 1) every row has exactly 11 tab-separated fields
awk -F'\t' 'NF!=11{print "bad fields on line " NR ": " NF; exit 1}' orchard_database_NE_states_021025_formatted.txt

# 2) no double-tabs (which would imply empty fields)
grep -nP '\t\t' orchard_database_NE_states_021025_formatted.txt | head


head -20 orchard_database_011025.txt
head -20 orchard_database_NE_states_021025_formatted.txt

# copy out to storage
cp orchard_database_NE_states_021025.txt /lustre2/home/lc736_0001/diet/songbird_coi_database/final_outputs/
cp orchard_database_NE_states_021025_formatted.txt /lustre2/home/lc736_0001/diet/songbird_coi_database/final_outputs/



##############################################################################
#*  UPDATE R SCRIPT (gbif_filter_states.R) HAVE NOT DONE THIS DO THIS AND TEST BEFORE NEXT RUN
##############################################################################
##### UPDATE R SCRIPT WITH THIS BELOW:
# before writing, turn "" -> NA
for (j in seq_along(dt_output)) {
  if (is.character(dt_output[[j]])) {
    set(dt_output, which(dt_output[[j]] == ""), j, NA_character_)
  }
}

# exact CRABS format
fwrite(
  dt_output,
  file = output_path,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  col.names = FALSE,
  showProgress = FALSE
)
