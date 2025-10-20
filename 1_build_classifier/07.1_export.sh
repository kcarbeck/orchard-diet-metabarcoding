# export databases for specific projects
# katherine carbeck
# 11 sep 2025

##############################################################################
#*        QIIME2-compatible formats: sequences and taxonomy
##############################################################################
crabs --export \
  --input /workdir/kcarbeck/final_outputs/orchard_database_NE_states_021025_formatted.txt \
  --output /workdir/kcarbeck/final_outputs/qiime_orchard_sequences_NE_states_141025.fasta \
  --export-format 'qiime-fasta' &
# |            Function | Export CRABS database to QIIME-FASTA format
# |         Import data | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:00:04
# |      Exporting data | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:00:02
# |            Results | Written 1605261 sequences to /workdir/kcarbeck/final_outputs/qiime_orchard_sequences_NE_states_141025.fasta out of 1605261 initial sequences (100.0%)

crabs --export \
  --input /workdir/kcarbeck/final_outputs/orchard_database_NE_states_021025_formatted.txt \
  --output /workdir/kcarbeck/final_outputs/qiime_orchard_taxonomy_NE_states_141025.txt \
  --export-format 'qiime-text'
# |            Function | Export CRABS database to QIIME-TEXT format
# |         Import data | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:00:04
# |      Exporting data | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:00:02
# |              Results | Written 1605261 sequences to /workdir/kcarbeck/final_outputs/qiime_orchard_taxonomy_NE_states_141025.txt out of 1605261 initial sequences (100.0%)

cp /workdir/kcarbeck/final_outputs/qiime_orchard_sequences_NE_states_141025.fasta /lustre2/home/lc736_0001/diet/songbird_coi_database/final_outputs/
cp /workdir/kcarbeck/final_outputs/qiime_orchard_taxonomy_NE_states_141025.txt /lustre2/home/lc736_0001/diet/songbird_coi_database/final_outputs/


