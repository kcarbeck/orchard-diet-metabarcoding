# filter the database
# katherine carbeck
# 10 sep 2025

# initial dereplication (species-level)
crabs --dereplicate \
  --input processed/anml_aligned_recovered.txt \
  --output processed/anml_aligned_recovered_dereplicated.txt \
  --dereplication-method 'unique_species'
#  Results | Written 3945345 unique sequences to processed/anml_aligned_recovered_dereplicated.txt out of 12072121 initial sequences (32.68%)

# quality filtering for COI metabarcoding
# removed --no-species-id
crabs --filter \
  --input processed/anml_aligned_recovered_dereplicated.txt \
  --output processed/anml_aligned_recovered_dereplicated_filtered.txt \
  --minimum-length 150 \
  --maximum-length 230 \
  --maximum-n 6 \
  --environmental \
  --rank-na 5
# Results | Written 3184691 filtered sequences to processed/anml_aligned_recovered_dereplicated_filtered.txt out of 3945345 initial sequences (80.72%)

