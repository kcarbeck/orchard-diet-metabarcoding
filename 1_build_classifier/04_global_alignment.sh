# retrieve amplicons without primer-binding regions
# katherine carbeck
# 10 sep 2025


# dereplicate
crabs --dereplicate --input processed/full_merged.clean.txt \
  --output processed/dereplicated.txt \
  --dereplication-method 'unique_species'

#filter
crabs --filter --input processed/dereplicated.txt \
  --output processed/dereplicated_filtered.txt \
  --environmental \
  --no-species-id \
  --rank-na 3


##############################################################################
#*             Recover additional sequences via global alignment
##############################################################################
# It is common practice to remove primer-binding regions from reference sequences when deposited in an online database. Therefore, when the reference sequence was generated using the same forward and/or reverse primer as searched for in the --in-silico-pcr function, the --in-silico-pcr function will have failed to recover the amplicon region of the reference sequence. To account for this possibility, CRABS has the option to run a Pairwise Global Alignment, implemented using VSEARCH v 2.16.0, to extract amplicon regions for which the reference sequence does not contain the full forward and reverse primer-binding regions. To accomplish this, the --pairwise-global-alignment function takes in the originally downloaded database file using the --input parameter. The database to be searched against is the output file from the --in-silico-pcr and can be specified using the --amplicons parameter. The output file can be specified using the --output parameter. The primer sequences, only used to calculate basepair length, can be set with the --forward and --reverse parameters. As the --pairwise-global-alignment function can take a long time to run for large databases, sequence length can be restricted to speed up the process using the --size-select parameter. Minimum percentage identity and query coverage can be specified using the --percent-identity and --coverage parameters, respectively. --percent-identity should be provided as a percentage value between 0 and 1 (e.g., 95% = 0.95), while --coverage should be provided as a percentage value between 0 and 100 (e.g., 95% = 95). By default, the --pairwise-global-alignment function is restricted to retain sequences where primer sequences are not fully present in the reference sequence (alignment starting or ending within the length of the forward or reverse primer). When the --all-start-positions parameter is provided, positive hits will be included when the alignment is found outside the range of the primer-binding regions (missed by --in-silico-pcr function due to too many mismatches in the primer-binding region). We do not recommend using the --all-start-positions, as it is very unlikely a barcode will be amplified using the specified primer set of the --in-silico-pcr function when more than 4 mismatches are present in the primer-binding regions.


# Recover additional sequences without complete primer regions
crabs --pairwise-global-alignment \
  --input processed/dereplicated_filtered.txt \
  --amplicons processed/anml_amplicons_primary.txt \
  --output processed/anml_aligned_recovered.txt \
  --forward GGTCAACAAATCATAAAGATATTGG \
  --reverse GGWACTAATCAATTTCCAAATCC \
  --percent-identity 0.92 \
  --coverage 85 \
  --size-select 2500 \
  --threads 20
#|            Function | Retrieve amplicons without primer-binding regions
#|         Import data | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:00:56
#|  Transform to fasta | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:00:27
#|  Pairwise alignment | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 11:38:20
#|Parse alignment data | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:00:00
#|      Exporting data | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:00:30
#|             Results | Retrieved 56199 amplicons without primer-binding regions from 737092 sequences

#plot the amplicon length distribution
crabs --amplicon-length-figure \
  --input processed/anml_aligned_recovered.txt \
  --output plots/amplicon_length.png \
  --tax-level 4
