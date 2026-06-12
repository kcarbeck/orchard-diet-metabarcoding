#!/bin/bash
# =============================================================================
# merge plates and years, then run decontam
# =============================================================================
# author: katherine carbeck
# date: nov 2025
#
# this script merges all sequencing plates/years into one dataset,
# then runs decontam to identify and remove contaminants
#
# pipeline order:
#   1. merge metadata (creates Year_Plate batch identifier)
#   2. merge feature tables from all plates
#   3. merge rep-seqs from all plates
#   4. merge taxonomy files
#   5. run decontam on merged data with batch mode
#   6. filter out contaminants
#
# why merge before decontam?
#   - more statistical power (more negatives = better contaminant detection)
#   - batch mode accounts for plate-to-plate variation
#   - simpler pipeline (one decontam run instead of multiple)
#
# =============================================================================

# =============================================================================
# configuration - edit these paths for your data
# =============================================================================

# working directory (where your qiime2 files are)
WORK_DIR="/lustre2/home/lc736_0001/orchards"
cd "$WORK_DIR"

# -----------------------------------------------------------------------------
# input files
# -----------------------------------------------------------------------------
# raw feature tables before decontam
TABLE_2024_P1="$WORK_DIR/2024_rerun3/table_plate1.qza"
TABLE_2025_P1="$WORK_DIR/2025/2025/table_plate1.qza"
TABLE_2025_P2="$WORK_DIR/2025_redo/table_plate2.qza"

# rep-seqs files 
REPSEQS_2024_P1="$WORK_DIR/2024_rerun3/rep-seqs_plate1.qza"
REPSEQS_2025_P1="$WORK_DIR/2025/2025/rep-seqs_plate1.qza"
REPSEQS_2025_P2="$WORK_DIR/2025_redo/rep-seqs_plate2.qza"

# taxonomy files (one per year, or one per plate - adjust as needed)
TAXONOMY_2024="$WORK_DIR/2024_rerun3/classifier_out/nb_classified_taxonomy_111125.qza"
TAXONOMY_2025_P1="$WORK_DIR/2025/classifier_out/nb_classified_taxonomy_141025.qza"
TAXONOMY_2025_P2="$WORK_DIR/2025_redo/classifier_out/nb_classified_taxonomy_112425.qza"

# metadata files
META_2024="$WORK_DIR/2024_rerun3/2024_metadata_3.txt"
META_2025_P1="$WORK_DIR/2025/2025/fecal_sampling_metadata_plate1.txt"
META_2025_P2="$WORK_DIR/2025_redo/fecal_sampling_metadata_plate2.txt"

# output directory (use absolute path to avoid cwd issues)
OUT_DIR="$WORK_DIR/merged"
mkdir -p "$OUT_DIR"

# ==============================================
# step 1: merge metadata"
# ==============================================
# this step combines all metadata files and creates:
#     - Year_Plate column for batch identification
#     - control_role column (sample, neg, pos)
#     - unique ids for duplicate controls
# 

# the r script 04.1_merge_metadata.R (already run) handles:
#   - normalizing column names (#SampleID, sample-id, etc.)
#   - creating Year_Plate column
#   - prefixing duplicate control ids
#   - generating rename maps if needed
 
# -----------------------------------------------------------------------------
# collect inputs from configured variables above (no external lists assumed)
# -----------------------------------------------------------------------------
TABLE_FILES=()
REPSEQS_FILES=()
TAXONOMY_FILES=()

echo "  collecting configured inputs..."
# feature tables
for f in \
    "$TABLE_2024_P1" \
    "$TABLE_2025_P1" \
    "$TABLE_2025_P2" \
; do
    [ -n "$f" ] && [ -f "$f" ] && TABLE_FILES+=("$f")
done
# rep-seqs
for f in \
    "$REPSEQS_2024_P1" \
    "$REPSEQS_2025_P1" \
    "$REPSEQS_2025_P2" \
; do
    [ -n "$f" ] && [ -f "$f" ] && REPSEQS_FILES+=("$f")
done
# taxonomy
for f in \
    "$TAXONOMY_2024" \
    "$TAXONOMY_2025_P1" \
    "$TAXONOMY_2025_P2" \
; do
    [ -n "$f" ] && [ -f "$f" ] && TAXONOMY_FILES+=("$f")
done
 
printf '%s\n' "${TABLE_FILES[@]}"
printf '%s\n' "${REPSEQS_FILES[@]}"
printf '%s\n' "${TAXONOMY_FILES[@]}"

# -----------------------------------------------------------------------------
# apply rename maps to feature tables 
# -----------------------------------------------------------------------------
# If files exist matching rename_maps/rename_map_*.tsv, apply ALL maps
# sequentially to EACH configured table. Only matching sample IDs
# are changed; others remain unchanged. Outputs *_renamed.qza into OUT_DIR.


#==============================================
# check for rename maps
# ==============================================
echo "  looking in: $OUT_DIR/rename_maps/"
echo "  pattern: rename_map_*.tsv"
echo "  TABLE_FILES has ${#TABLE_FILES[@]} entries"

if [ -d "$OUT_DIR/rename_maps" ]; then
    echo "  directory exists: YES"
    echo "  files found:"
    ls -la "$OUT_DIR/rename_maps/" 2>/dev/null || echo "    (none or cannot list)"
else
    echo "  directory exists: NO"
    echo "  skipping rename step (create $OUT_DIR/rename_maps/ with rename_map_*.tsv files if needed)"
fi

if compgen -G "$OUT_DIR/rename_maps/rename_map_*.tsv" > /dev/null 2>&1; then
    echo "  rename maps detected in $OUT_DIR/rename_maps; applying to feature tables..."
    RENAMED_TABLE_FILES=()
    mkdir -p "$OUT_DIR"
    for table_path in "${TABLE_FILES[@]}"; do
        # Include parent directory in output name to avoid overwrites
        # e.g., 2024_rerun3/table_plate1.qza -> 2024_rerun3_table_plate1_renamed.qza
        parent_dir="$(basename "$(dirname "$table_path")")"
        base="$(basename "$table_path" .qza)"
        output_name="${parent_dir}_${base}"
        
        # Determine which year this table belongs to based on path
        if [[ "$table_path" == *"2024"* ]]; then
            year_pattern="2024"
        elif [[ "$table_path" == *"2025"* ]]; then
            year_pattern="2025"
        else
            year_pattern=""  # apply all maps if year not detected
        fi
        
        current_input="$table_path"
        temp_dir="$(mktemp -d 2>/dev/null || mktemp -d -t tmp)"
        maps_applied=0
 
        for map_file in "$OUT_DIR"/rename_maps/rename_map_*.tsv; do
            map_name="$(basename "$map_file")"
            
            # Only apply rename maps that match this table's year
            if [ -n "$year_pattern" ] && [[ "$map_name" != *"$year_pattern"* ]]; then
                continue  # skip maps from different years
            fi
            
            next_output="$temp_dir/${output_name}_tmp.qza"
            qiime feature-table rename-ids \
                --i-table "$current_input" \
                --m-metadata-file "$map_file" \
                --m-metadata-column new_id \
                --p-axis sample \
                --o-renamed-table "$next_output"
            current_input="$next_output"
            ((maps_applied++))
            echo "      applied: $map_name"
        done
 
        final_output="$OUT_DIR/${output_name}_renamed.qza"
        cp -f "$current_input" "$final_output"
        rm -rf "$temp_dir"
        RENAMED_TABLE_FILES+=("$final_output")
        echo "    wrote: $final_output ($maps_applied maps applied)"
    done
 
    # Use renamed tables for downstream merge
    if [ "${#RENAMED_TABLE_FILES[@]}" -gt 0 ]; then
        TABLE_FILES=("${RENAMED_TABLE_FILES[@]}")
    fi
else
    echo "  no rename maps found matching pattern; skipping rename step"
fi


# ==============================================
# step 2: merge feature tables"
# ==============================================
# combine asv count tables from all plates

#the metadata file has crazy headers so fix that first:
# backup metadata file
cp "$OUT_DIR/all_plates_metadata.tsv" "$OUT_DIR/all_plates_metadata.tsv.bak"
# replace the header line with the correct one (all 17 columns)
sed -i '1s/.*/SampleID\tcontrol_role\tYear\tPlateName\tYear_Plate\tWell\tYear_orig\tSpecies\tSite\tDate\tAge\tSex\tSubstrate_group\tPoop_texture\tPoop_color\tLetter_codes\tBand_Number/' "$OUT_DIR/all_plates_metadata.tsv"
#check header
head -2 all_plates_metadata.tsv | cat -A   # shows tabs as ^I
awk -F'\t' '{print NF}' all_plates_metadata.tsv | sort | uniq -c  # should show all rows have 17 cols
 
if [ "${#TABLE_FILES[@]}" -eq 0 ]; then
    echo "ERROR: No feature table files found. Provide paths in $TABLES_LIST" >&2
    exit 1
fi
 
if [ "${#TABLE_FILES[@]}" -eq 1 ]; then
    echo "  only one table provided; copying to $OUT_DIR/table_merged.qza"
    cp -f "${TABLE_FILES[0]}" "$OUT_DIR/table_merged.qza"
else
    MERGE_TABLE_ARGS=()
    for f in "${TABLE_FILES[@]}"; do
        MERGE_TABLE_ARGS+=( --i-tables "$f" )
    done
    qiime feature-table merge \
        "${MERGE_TABLE_ARGS[@]}" \
        --o-merged-table "$OUT_DIR/table_merged.qza"
fi

echo "  wrote: $OUT_DIR/table_merged.qza"

# create summary visualization
qiime feature-table summarize \
    --i-table "$OUT_DIR/table_merged.qza" \
    --m-sample-metadata-file all_plates_metadata.tsv \
    --o-visualization "$OUT_DIR/table_merged.qzv"


#==============================================
# step 3: merge rep-seqs
#==============================================
# combine representative sequences from all plates

 
if [ "${#REPSEQS_FILES[@]}" -eq 0 ]; then
    echo "ERROR: No rep-seqs files found. Provide paths in $REPSEQS_LIST" >&2
    exit 1
fi
 
if [ "${#REPSEQS_FILES[@]}" -eq 1 ]; then
    echo "  only one rep-seqs provided; copying to $OUT_DIR/rep-seqs_merged.qza"
    cp -f "${REPSEQS_FILES[0]}" "$OUT_DIR/rep-seqs_merged.qza"
else
    MERGE_REPSEQS_ARGS=()
    for f in "${REPSEQS_FILES[@]}"; do
        MERGE_REPSEQS_ARGS+=( --i-data "$f" )
    done
    qiime feature-table merge-seqs \
        "${MERGE_REPSEQS_ARGS[@]}" \
        --o-merged-data "$OUT_DIR/rep-seqs_merged.qza"
fi

echo "  wrote: $OUT_DIR/rep-seqs_merged.qza"

qiime tools export \
    --input-path rep-seqs_merged.qza \
    --output-path exported_seqs

#==============================================
# step 4: merge taxonomy
#==============================================
# combine taxonomy annotations from all plates

# export taxonomy files to tsv format
if [ "${#TAXONOMY_FILES[@]}" -eq 0 ]; then
    echo "WARNING: No taxonomy files found; skipping taxonomy export/merge." >&2
else
    for t in "${TAXONOMY_FILES[@]}"; do
        base="$(basename "$t" .qza)"
        out_dir="$OUT_DIR/export_${base}"
        mkdir -p "$out_dir"
        qiime tools export --input-path "$t" --output-path "$out_dir"
        echo "  exported: $t -> $out_dir"
    done
fi

if [ "${#TAXONOMY_FILES[@]}" -gt 0 ]; then
    # run r script to merge taxonomy
    # this keeps the deepest classification rank and highest confidence
    R_MERGE_TAX="helper_scripts/merge_taxonomy.R"
    if [ -f "$R_MERGE_TAX" ]; then
        Rscript "$R_MERGE_TAX"
    else
        echo "WARNING: $R_MERGE_TAX not found; skipping taxonomy merge." >&2
    fi
 
    # import merged taxonomy back to qiime2
    qiime tools import \
        --type 'FeatureData[Taxonomy]' \
        --input-path taxonomy_merged.tsv \
        --output-path "$OUT_DIR/taxonomy_merged.qza"
 
    echo "  wrote: $OUT_DIR/taxonomy_merged.qza"
else
    echo "  skipped taxonomy merge/import (no taxonomy inputs found)"
fi
#total rows before merging: 43511
#unique features: 40294
#rows after merging: 40294


#==============================================
# step 5: validate merged metadata
#==============================================

qiime metadata tabulate \
    --m-input-file all_plates_metadata.tsv \
    --o-visualization "$OUT_DIR/all_plates_metadata.qzv"

# open this in qiime2 view to check your metadata, looks good!


#==============================================
# merge complete!
#==============================================

# what to do next:
#  1. run decontam: 05_decontam_prevalence.R
#  2. review threshold_sweep.csv and diagnostic plots
#  3. run the filter commands shown above
#  4. proceed to 06_qc_and_get_analysis_files.sh

# output files in $OUT_DIR:
#  - table_merged.qza         (combined feature table)
#  - rep-seqs_merged.qza      (combined rep sequences)
#  - taxonomy_merged.qza      (combined taxonomy)
#  - all_plates_metadata.qzv  (metadata visualization)

# ==============================================
