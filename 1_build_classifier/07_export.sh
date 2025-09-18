# export databases for specific projects
# katherine carbeck
# 11 sep 2025

##############################################################################
#*        QIIME2-compatible formats: sequences and taxonomy
##############################################################################
crabs --export \
  --input final_outputs/orchard_database.txt \
  --output final_outputs/qiime_orchard_sequences.fasta \
  --export-format 'qiime-fasta'

crabs --export \
  --input final_outputs/orchard_database.txt \
  --output final_outputs/qiime_orchard_taxonomy.txt \
  --export-format 'qiime-text'

##############################################################################
#*        Additional classifier formats: sintax and blast-tax
##############################################################################
crabs --export \
  --input final_outputs/orchard_database.txt \
  --output final_outputs/sintax_orchard_database.fasta \
  --export-format sintax

crabs --export \
  --input final_outputs/orchard_database.txt \
  --output final_outputs/blast_orchard_database \
  --export-format blast-tax
