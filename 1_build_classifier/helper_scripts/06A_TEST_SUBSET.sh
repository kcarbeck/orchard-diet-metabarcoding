

##############################################################################
#*                   TEST DATABASE
##############################################################################
# initial dereplication (species-level)
crabs --dereplicate \
  --input processed/anml_amplicons_primary.txt \
  --output test/test_dereplicated.txt \
  --dereplication-method 'unique_species'

# quality filtering for COI metabarcoding
crabs --filter \
  --input test/test_dereplicated.txt \
  --output test/test_dereplicated_filtered.txt \
  --minimum-length 150 \
  --maximum-length 210 \
  --maximum-n 6 \
  --environmental \
  --no-species-id \
  --rank-na 3

 crabs --subset \
 --input test/test_dereplicated_filtered.txt \
 --include 'Lepidoptera;Diptera;Coleoptera;Hymenoptera;Hemiptera;Araneae;Collembola;Orthoptera;Acari;Neuroptera;Dermaptera;Thysanoptera' \
 --output test/test_database.txt

crabs --diversity-figure --input test/test_database.txt --output test/test_diversity_4.png --tax-level 4

crabs --diversity-figure --input test/test_database.txt --output test/test_diversity_3.png --tax-level 3

crabs --diversity-figure --input test/test_database.txt --output test/test_diversity_2.png --tax-level 2

#plot the amplicon length distribution
crabs --amplicon-length-figure \
  --input test/test_database.txt \
  --output test/test_amplicon_length.png \
  --tax-level 4
###!
crabs --phylogenetic-tree --input test/test_database.txt --output test/test_phylo --tax-level 4 --species 'Carcharodon carcharias+Squalus acanthias'