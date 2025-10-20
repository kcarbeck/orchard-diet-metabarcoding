# dereplicate reference database
# katherine carbeck
# 14 oct 2025

# import sequences into qiime
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /workdir/kcarbeck/final_outputs/qiime_orchard_sequences_NE_states_141025.fasta \
  --output-path /workdir/kcarbeck/final_outputs/qiime_orchard_sequences_NE_states_141025.qza 

# import taxonomy
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /workdir/kcarbeck/final_outputs/qiime_orchard_taxonomy_NE_states_141025.txt \
  --output-path /workdir/kcarbeck/final_outputs/qiime_orchard_taxonomy_NE_states_141025.qza 


cp /workdir/kcarbeck/final_outputs/qiime_orchard_sequences_NE_states_141025.qza /lustre2/home/lc736_0001/diet/songbird_coi_database/final_outputs/
cp /workdir/kcarbeck/final_outputs/qiime_orchard_taxonomy_NE_states_141025.qza /lustre2/home/lc736_0001/diet/songbird_coi_database/final_outputs/

##############################################################################
#*               2. dereplicate using 2 methodss:
#*                       super mode & LCA 
##############################################################################
# "super" finds the LCA consensus while giving preference to majority labels and collapsing substrings into superstrings. For example, when a more specific taxonomy does not contradict a less specific taxonomy, the more specific is chosen. That is, "g__Faecalibacterium; s__prausnitzii", will be preferred over "g__Faecalibacterium; s__" 

qiime rescript dereplicate \
  --i-sequences /workdir/kcarbeck/final_outputs/qiime_orchard_sequences_NE_states_141025.qza \
  --i-taxa /workdir/kcarbeck/final_outputs/qiime_orchard_taxonomy_NE_states_141025.qza \
  --p-mode super \
  --p-threads 22 \
  --o-dereplicated-sequences /workdir/kcarbeck/final_outputs/ref_derep_seqs_141025.qza \
  --o-dereplicated-taxa /workdir/kcarbeck/final_outputs/ref_derep_tax_141025.qza

qiime rescript dereplicate \
  --i-sequences /workdir/kcarbeck/final_outputs/qiime_orchard_sequences_NE_states_141025.qza \
  --i-taxa /workdir/kcarbeck/final_outputs/qiime_orchard_taxonomy_NE_states_141025.qza \
  --p-mode lca \
  --p-threads 22 \
  --o-dereplicated-sequences /workdir/kcarbeck/final_outputs/ref_derep_lca_seqs_141025.qza \
  --o-dereplicated-taxa /workdir/kcarbeck/final_outputs/ref_derep_lca_tax_141025.qza

cp /workdir/kcarbeck/final_outputs/ref_derep_seqs_141025.qza /lustre2/home/lc736_0001/diet/songbird_coi_database/final_outputs/
cp /workdir/kcarbeck/final_outputs/ref_derep_tax_141025.qza /lustre2/home/lc736_0001/diet/songbird_coi_database/final_outputs/
cp /workdir/kcarbeck/final_outputs/ref_derep_lca_seqs_141025.qza /lustre2/home/lc736_0001/diet/songbird_coi_database/final_outputs/
cp /workdir/kcarbeck/final_outputs/ref_derep_lca_tax_141025.qza /lustre2/home/lc736_0001/diet/songbird_coi_database/final_outputs/


##############################################################################
#*               2. export sequences and taxonomy to text files
##############################################################################
qiime tools peek /workdir/kcarbeck/final_outputs/ref_derep_seqs_141025.qza
# expect: type: FeatureData[Sequence]

outdir=/workdir/kcarbeck/final_outputs
mkdir -p "$outdir"

qiime tools export \
  --input-path /workdir/kcarbeck/final_outputs/ref_derep_seqs_141025.qza \
  --output-path "$outdir"

# export to fasta
mv "$outdir"/dna-sequences.fasta /workdir/kcarbeck/final_outputs/ref_derep_seqs_141025.fasta


# make a 2 column tsv with feature-id and sequence
awk '
  BEGIN{h=""; s=""}
  /^>/{
    if(h!=""){print h "\t" s}
    h=substr($0,2); s=""
    next
  }
  {s=s $0}
  END{if(h!=""){print h "\t" s}}
' /workdir/kcarbeck/final_outputs/ref_derep_seqs_141025.fasta \
  > /workdir/kcarbeck/final_outputs/ref_derep_seqs_141025.tsv


# taxonomy file:
qiime tools export \
  --input-path /workdir/kcarbeck/final_outputs/ref_derep_lca_tax_141025.qza \
  --output-path "$outdir"



