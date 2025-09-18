# subset databases for specific projects
# katherine carbeck
# 10 sep 2025

# filter database for relevant taxa
crabs --subset \
 --input processed/anml_aligned_recovered_dereplicated_filtered.txt \
 --include 'Lepidoptera;Diptera;Coleoptera;Hymenoptera;Hemiptera;Araneae;Collembola;Orthoptera;Acari;Neuroptera;Dermaptera;Thysanoptera' \
 --output final_outputs/orchard_database.txt
#  Results | Written 2845633 subsetted sequences to final_outputs/orchard_database.txt out of 3184691 initial sequences (89.35%)






