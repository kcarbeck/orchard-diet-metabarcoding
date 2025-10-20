# subset databases for specific projects
# katherine carbeck
# 10 sep 2025

# filter database for relevant taxa
crabs --subset \
 --input processed/anml_aligned_recovered_dereplicated_filtered_011025.txt \
 --include 'Lepidoptera;Diptera;Coleoptera;Hymenoptera;Hemiptera;Araneae;Collembola;Orthoptera;Acari;Neuroptera;Dermaptera;Thysanoptera' \
 --output final_outputs/orchard_database_011025.txt
# |             Results | Written 2677174 subsetted sequences to final_outputs/orchard_database_011025.txt out of 3023466 initial sequences (88.55%)

