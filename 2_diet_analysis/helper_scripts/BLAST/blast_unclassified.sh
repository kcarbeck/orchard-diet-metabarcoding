#!/bin/bash
# =============================================================================
# BLAST poorly classified sequences
# =============================================================================
# This script identifies ASVs that are only classified to Arthropoda/Insecta
# and BLASTs them against NCBI to see what they actually are
# =============================================================================

# activate qiime2
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source /programs/miniconda3/bin/activate qiime2-amplicon-2024.10

WORK_DIR="/lustre2/home/lc736_0001/orchards/merged"
cd "$WORK_DIR"

OUTPUT_DIR="$WORK_DIR/blast_unclassified"
mkdir -p "$OUTPUT_DIR"

# =============================================================================
# STEP 1: Export taxonomy and rep-seqs
# =============================================================================

# export taxonomy to TSV
qiime tools export \
    --input-path taxonomy_merged.qza \
    --output-path "$OUTPUT_DIR/taxonomy_export"

# export representative sequences to FASTA
qiime tools export \
    --input-path rep-seqs_merged.qza \
    --output-path "$OUTPUT_DIR/seqs_export"

# export feature table to get read counts
qiime tools export \
    --input-path analysis_tables/table_arthropoda.qza \
    --output-path "$OUTPUT_DIR/table_export"

# convert biom to tsv
biom convert \
    -i "$OUTPUT_DIR/table_export/feature-table.biom" \
    -o "$OUTPUT_DIR/feature_table.tsv" \
    --to-tsv


# =============================================================================
# STEP 2: Identify poorly classified ASVs (run in R or manually)
# =============================================================================
# The taxonomy file is at: $OUTPUT_DIR/taxonomy_export/taxonomy.tsv
# Format: Feature ID, Taxon, Confidence
#
# Look for ASVs where Taxon ends at Arthropoda or Insecta level
# i.e., lots of NA or "Unassigned" at lower ranks
#

# GO TO R SCRIPT: identify_blast_unclassified.R
# then come back here to blast

# =============================================================================
# STEP 3: Identify poorly classified ASVs and prepare for BLAST
# =============================================================================
# use output from r to blast 
# output is in: /lustre2/home/lc736_0001/orchards/merged/blast_unclassified

# FIRST, set up working directory with local nt database
mkdir -p /workdir/$USER
cd /workdir/$USER

# copy nt database (this might take a few minutes - it's large)
cp /shared_data/genome_db/BLAST_NCBI/nt* ./

# copy taxdb for species names (optional but helpful)
cp /shared_data/genome_db/BLAST_NCBI/taxdb* ./ 2>/dev/null || echo "taxdb not available"

# copy your query file
cp /lustre2/home/lc736_0001/orchards/merged/blast_unclassified/all_poorly_classified.fasta ./

# run local blast (with multiple threads - it will be much faster than remote blast)
blastn -query all_poorly_classified.fasta \
       -db nt \
       -num_threads 16 \
       -max_target_seqs 5 \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames' \
       -out blast_all_results.tsv


# How many sequences got hits?
wc -l blast_all_results.tsv
# 49777 blast_all_results.tsv

# Look at top hits
head -20 blast_all_results.tsv

# Count unique species found
cut -f14 blast_all_results.tsv | sort | uniq -c | sort -rn | head -30
#   1706 Rhagio tringarius
#   1047 Barypeithes pellucidus
#   1020 Trachelipus rathkii
#    993 Amphipoea velata
#    955 Sylvicola alternatus
#    685 Philoscia muscorum
#    671 Macaca fascicularis
#    528 Severe acute respiratory syndrome coronavirus 2
#    496 Ichneumoninae sp. BOLD-2016
#    441 Avena sativa
#    425 Porcellio spinicornis
#    399 uncultured Arcellinida
#    342 Phaeocystis globosa
#    335 Valenzuela flavidus
#    299 Auxenochlorella protothecoides
#    299 Amphipsocidae sp. BOLD:AAH3229
#    286 Carex neurocarpa
#    284 Cylindroiulus punctatus
#    278 Carex depauperata
#    270 Adineta vaga
#    260 Pleuretra lineata
#    251 Zanclognatha pedipilalis
#    246 Carex divulsa
#    239 Chaetocnema concinna
#    230 Ixodes scapularis
#    229 Homo sapiens
#    229 Cecidomyiidae sp. BOLD-2016
#    220 Acanthamoeba sp.
#    209 Bangia atropurpurea
#    204 Peloptulus phaenotus

###############
# AFTER DONE:
# copy out to storage:
cp blast_all_results.tsv /lustre2/home/lc736_0001/orchards/merged/blast_unclassified/

# clean up 
rm -rf /workdir/kcarbeck/nt*


# =============================================================================
# STEP 4: Merge BLAST results with taxonomy
# =============================================================================
# use merge_blast_with_taxonomy.R to merge blast results with taxonomy
# output is in: /lustre2/home/lc736_0001/orchards/merged/blast_unclassified/enhanced_taxonomy

# DONE WIT THIS TANGENT --> NOW GO BACK TO 06_qc_and_get_analysis_files.sh STEP 8