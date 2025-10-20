# filter the database
# katherine carbeck
# 10 sep 2025

# initial dereplication (species-level)
crabs --dereplicate \
  --input processed/anml_aligned_recovered.txt \
  --output processed/anml_aligned_recovered_dereplicated.txt \
  --dereplication-method 'unique_species'
#  Results | Written 3945345 unique sequences to processed/anml_aligned_recovered_dereplicated.txt out of 12072121 initial sequences (32.68%)

# quality filtering for COI 
# removed --no-species-id
crabs --filter \
  --input processed/anml_aligned_recovered_dereplicated.txt \
  --output processed/anml_aligned_recovered_dereplicated_filtered_011025.txt \
  --minimum-length 160 \
  --maximum-length 210 \
  --maximum-n 2 
#|             Results | Written 3023466 filtered sequences to
#processed/anml_aligned_recovered_dereplicated_filtered_011025.txt out of 3945345 initial sequences (76.63%)
#|                     | Minimum length filter: 371275 sequences not passing filter (9.41%)
#|                     | Maximum length filter: 111230 sequences not passing filter (2.82%)
#|                     | Maximum ambiguous bases filter: 504905 sequences not passing filter (12.8%)