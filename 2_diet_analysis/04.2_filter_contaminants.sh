# filter contaminants from per-plate/year tables
# katherine carbeck
# 20 oct 2025

# now that we've identified the contaminants using the negatives, we drop those ASVs
qiime feature-table filter-features \
  --i-table /workdir/kcarbeck/2025/table_plate1.qza \
  --m-metadata-file /workdir/kcarbeck/2025/decontam_plate1/contaminant_feature_ids_thr_010.txt \
  --p-exclude-ids \
  --o-filtered-table /workdir/kcarbeck/2025/table_plate1_nocontam.qza

# summarize the cleaned table and look at read counts per sample
qiime feature-table summarize \
  --i-table /workdir/kcarbeck/2025/table_plate1_nocontam.qza \
  --o-visualization /workdir/kcarbeck/2025/table_plate1_nocontam.qzv

# filter rep seqs to keep only the seqs still present in the filtered table so the IDs match for later
qiime feature-table filter-seqs \
  --i-data /workdir/kcarbeck/2025/rep-seqs_plate1.qza \
  --i-table /workdir/kcarbeck/2025/table_plate1_nocontam.qza \
  --o-filtered-data /workdir/kcarbeck/2025/rep-seqs_plate1_nocontam.qza
# you should still have negatives and mocks in the table at this step because we want to keep them for provenance/QC and we'll remove them after merging all the plates/batches

# Saved FeatureTable[Frequency] to: /workdir/kcarbeck/2025/table_plate1_nocontam.qza
# Saved Visualization to: /workdir/kcarbeck/2025/table_plate1_nocontam.qzv
#Saved FeatureData[Sequence] to: /workdir/kcarbeck/2025/rep-seqs_plate1_nocontam.qza

# copied files out to storage for tasos
cp -r decontam_plate1 /lustre2/home/lc736_0001/orchards/2025/
cp table_plate1_nocontam* /lustre2/home/lc736_0001/orchards/2025/decontam_plate1
cp rep-seqs_plate1_nocontam.qza /lustre2/home/lc736_0001/orchards/2025/decontam_plate1



# set the correct group everywhere
chgrp -R lc736_0001 /lustre2/home/lc736_0001/orchards/2025/decontam_plate1

# directories: 2770 (rwx for owner/group, setgid, no others)
find /lustre2/home/lc736_0001/orchards/2025/decontam_plate1 -type d -exec chmod 2770 {} +

# files: 770 (rwx for owner/group, no others)
find /lustre2/home/lc736_0001/orchards/2025/decontam_plate1 -type f -exec chmod 0770 {} +