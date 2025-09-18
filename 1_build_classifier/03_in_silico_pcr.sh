# extract amplicon regions through in silico PCR
# katherine carbeck
# 09 sep 2025

# CRABS extracts the amplicon region of the primer set by conducting an in silico PCR using the --in-silico-pcr function. CRABS uses cutadapt v 4.4 for the in silico PCR to increase speed of execution of traditional python code. Both the forward and reverse primer should be provided in 5'-3' direction using the '--forward' and '--reverse' parameters, respectively. CRABS will reverse complement the reverse primer. From version v 1.0.0, CRABS is capable of retaining barcodes in both direction using a single in silico PCR analysis. Hence, no reverse complementation step and rerunning of the in silico PCR is conducted, thereby significantly increasing execution speed. To retain sequences for which no primer-binding regions could be found, an output file can be specified for the --untrimmed parameter. The maximum allowed number of mismatches found in the primer-binding regions can be specified using the --mismatch parameter, with a default setting of 4. Finally, the in silico PCR analysis can be multithreaded in CRABS. By default the maximum number of threads are being used, but users can specify the number of threads to use with the --threads parameter.


conda activate CRABS

# first, make sure cutadapt is the version CRABS expects
cutadapt --version
# not 4.4, pin it (and avoid xopen>=2 which drags Windows-only metadata)
conda install -c bioconda "cutadapt=4.4" "xopen<2"

#back up and clean input (strip the BOLD HTML and any non-DNA sequences)
cp processed/full_merged.txt processed/full_merged.orig.txt

# drop rows containing obvious HTML/PHP error text
grep -Ev 'Fatal error|HTTP_Request2_MessageException|<html|<table|KOHANA|v3\.boldsystems\.org' \
  processed/full_merged.orig.txt > processed/full_merged.nohtml.txt

# keep only rows where the LAST column looks like a DNA sequence (IUPAC DNA letters)
# (Assumes the last column is the sequence in your CRABS-format TSV)
awk -F'\t' 'BEGIN{OFS="\t"} NR==1 || $NF ~ /^[ACGTRYSWKMBDHVN\-]+$/ {print}' \
  processed/full_merged.nohtml.txt > processed/full_merged.clean.txt

# see which IDs got dropped
comm -13 <(cut -f1 processed/full_merged.clean.txt | sort) \
          <(cut -f1 processed/full_merged.orig.txt  | sort)


##############################################################################
#*                 in silico PCR
##############################################################################
# Primary in silico PCR with relaxed mode (critical for COI recovery)
    # --mismatch 4: Optimal balance for arthropod COI diversity
    # ---relaxed can be specified. When used, CRABS will rerun the in silico PCR analysis after the first attempt, but now by only checking the presence of a single primer-binding region. If either the forward or reverse primer-binding region is found in this second attempt, the amplicon will be added to the --output file. This setting can be useful for amplicons where either the forward or reverse primer is commonly used as a barcoding primer.
    # --untrimmed: Captures sequences for alignment-based recovery
crabs --in-silico-pcr \
  --input processed/full_merged.clean.txt \
  --output processed/anml_amplicons_primary.txt \
  --forward GGTCAACAAATCATAAAGATATTGG \
  --reverse GGWACTAATCAATTTCCAAATCC \
  --mismatch 4 \
  --relaxed \
  --untrimmed processed/anml_failed_primary.txt \
  --threads 20
# 
# |            Function | Extract amplicons through in silico PCR
# |         Import data | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:01:14
# |  Transform to fasta | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:00:42
# |       In silico PCR | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:03:34
# | Transform untrimmed | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:01:57
# |      relaxed IS PCR | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:04:22
# |      Exporting data | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:00:27
# |      Exporting data | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00 0:00:41
# |             Results | Extracted 12015922 amplicons from 13851704 sequences (86.75%)
# |             Results | 11933264 amplicons were extracted by only the forward or reverse primer (99.31%)
