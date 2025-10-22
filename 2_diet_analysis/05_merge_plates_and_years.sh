# merge cleaned tables within each year
# katherine carbeck
# 20 oct 2025

########################################################
#*  EDIT FILE NAMES AND PATHS 
#* you will need to ensure that the sample names 
#* are unique across plates/years --- very important!!
########################################################


########################################################
#* export per-plate taxonomy to tsv
########################################################
## plate tags you expect to process (used only if rename maps exist)
PLATE_TAGS="2024A 2024B 2025A 2025B"


for t in taxonomy_2024A.qza taxonomy_2024B.qza taxonomy_2025A.qza taxonomy_2025B.qza; do
  out="export_${t%.qza}"
  qiime tools export --input-path "$t" --output-path "$out"
done
# you'll now have export_*/*taxonomy.tsv files


########################################################
#* R script to merge taxonomy rows 
#* (keep deepest rank, then highest confidence on ties)
########################################################
# look at the rscript before running to make sure file names/paths are correct
# may need to update some things within
Rscript 2_diet_analysis/merge_taxonomy.R



########################################################
#* import merged taxonomy back to qiime2
########################################################
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path taxonomy_merged.tsv \
  --output-path taxonomy_merged.qza

########################################################
#* build merged metadata file 
########################################################
# plate tags to use in IDs: 2024A, 2024B, 2025A, 2025B
#example per-plate metadata (qiime2-style TSV):
#  meta_2024A.tsv
#  meta_2024B.tsv
#  meta_2025A.tsv
#  meta_2025B.tsv

# merge metadata + enforce uniqueness + control_role (the control_role column will be useful later for filtering)

# note: there is a 'qiime metadata merge' command that can be used to merge 2 metadata files at a time, but fails when you have overlapping IDs and overlapping columns with conflicting values
# note2: DO NOT rename feature IDs since it will break links with taxonomy and rep-seqs

# look at the rscript before running to make sure file names/paths are correct
# may need to update some things within
Rscript 2_diet_analysis/merge_metadata.R

########################################################
#* important: if duplicates were found…
# You’ll see rename_maps/rename_map_<plate>.tsv files. Before merging any QIIME artifacts, rename the sample IDs in each per-plate table so they match the merged metadata
# if no duplicates were found, skip this step
########################################################

# run only if the R script warned about duplicates and produced rename_maps/*
if [ -d rename_maps ] && compgen -G "rename_maps/rename_map_*.tsv" > /dev/null; then
  for tag in ${PLATE_TAGS}; do
    if [ -f "rename_maps/rename_map_${tag}.tsv" ]; then
      qiime feature-table rename-ids \
        --i-table table_${tag}.qza \
        --m-metadata-file rename_maps/rename_map_${tag}.tsv \
        --m-metadata-column new_id \
        --p-axis sample \
        --o-renamed-table table_${tag}_renamed.qza
      # keep rep-seqs aligned with the renamed table
      qiime feature-table filter-seqs \
        --i-data  rep-seqs_${tag}.qza \
        --i-table table_${tag}_renamed.qza \
        --o-filtered-data rep-seqs_${tag}_renamed.qza
    fi
  done
fi

########################################################
#*  merge per plate tables & seqs
########################################################

qiime feature-table merge \
  --i-tables table_2024_plate1_nocontam.qza \
              table_2024_plate2_nocontam.qza \
              table_2025_plate1_nocontam.qza \
              table_2025_plate2_nocontam.qza \
  --o-merged-table table_merged_nocontam.qza

qiime feature-table merge-seqs \
  --i-data   rep-seqs_2024_plate1_nocontam.qza \
             rep-seqs_2024_plate2_nocontam.qza \
             rep-seqs_2025_plate1_nocontam.qza \
             rep-seqs_2025_plate2_nocontam.qza \
  --o-merged-data rep-seqs_merged_nocontam.qza


########################################################
#* validate merged metadata
########################################################
qiime metadata tabulate \
  --m-input-file all_plates_metadata.tsv \
  --o-visualization all_plates_metadata.qzv

########################################################
#* check coverage: do we have taxonomy for every merged feature?
########################################################
# list feature ids in merged rep-seqs
qiime tools export --input-path rep-seqs_merged_nocontam.qza --output-path export_repseqs
grep '^>' export_repseqs/dna-sequences.fasta | sed 's/^>//' | sort > ids_repseqs.txt

# list taxonomy feature ids
cut -f1 taxonomy_merged.tsv | tail -n +2 | sort > ids_taxonomy.txt

comm -23 ids_repseqs.txt ids_taxonomy.txt > missing_tax_ids.txt  # features lacking taxonomy
wc -l missing_tax_ids.txt
