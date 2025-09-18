# train QIIME2 classifier
# katherine carbeck
# 11 sep 2025

# activate qiime2 environment
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source /programs/miniconda3/bin/activate qiime2-amplicon-2024.10


##############################################################################
#*               1. import sequences and taxonomy
##############################################################################
# import sequences into qiime
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path final_outputs/qiime_orchard_sequences.fasta \
  --output-path final_outputs/qiime_orchard_sequences.qza &

# import taxonomy
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path final_outputs/qiime_orchard_taxonomy.txt \
  --output-path final_outputs/qiime_orchard_taxonomy.qza &


# how many refs am i training on?
grep -c '^>' final_outputs/qiime_orchard_sequences.fasta

# validate inputs
qiime tools validate final_outputs/qiime_orchard_sequences.qza
qiime tools validate final_outputs/qiime_orchard_taxonomy.qza


##############################################################################
#*                2.  VSEARCH CONSENSUS CLASSIFIER (PASS 1)
##############################################################################
# VSEARCH consensus classifier uses your FeatureData[Sequence] + FeatureData[Taxonomy] directly and is much lighter on RAM
# DOCUMENTATION: https://docs.qiime2.org/2024.10/plugins/available/feature-classifier/classify-consensus-vsearch/


# INPUT FILES:
# - orchards/rep-seqs_plate1.qza : your query ASVs (from DADA2/denoising)
# - final_outputs/qiime_orchard_sequences.qza : reference sequences (imported above)
# - final_outputs/qiime_orchard_taxonomy.qza : reference taxonomy (imported above)
#
# OUTPUT FILES (written into $OUTDIR):
# - classification.qza : per-ASV taxonomy assignments
# - search_results.qza : VSEARCH hit table used to make the assignments
#
# HOW TO RUN THIS STEP:
# 1) Activate QIIME2 (done at top of this script).
# 2) Verify the three INPUT files exist. If your query file is different (e.g., a
#    different plate or run), change --i-query to the correct .qza file.
# 3) Adjust the thresholds below (quick notes):
#      perc-identity (0.97): minimum % identity to count a hit in voting. higher favors precision (species-level), may reduce recall
#      weak-id (0.94): include weaker hits for higher-rank consensus only. relax to recover more genus/family.
#      query-cov (0.85): minimum fraction of the query covered by the hit. keep ≥ 0.8 for short amplicons to avoid partial/off-target matches.
#      maxaccepts (12): number of top hits to consider; more = stronger consensus.
#      maxrejects (256): how many non-matching hits to skip while searching.
#      min-consensus (0.8): fraction of considered hits that must agree per rank.
#      threads (22): CPUs to use; lower on laptops (e.g., 4–8).
# 4) Run this script. See the CHECK RESULTS section below to export/visualize.

# NOTES:
# if you want more species-level: try --p-perc-identity 0.985–0.99 and keep --p-min-consensus 0.8–0.9.
# if too many reads stop at high rank: lower --p-perc-identity slightly or increase --p-maxaccepts to 20.


# make a fresh, timestamped output dir
OUTDIR="orchards/vsearch_out_$(date +%Y%m%d_%H%M%S)"
# Train/classify:
qiime feature-classifier classify-consensus-vsearch \
  --i-query orchards/rep-seqs_plate1.qza \
  --i-reference-reads final_outputs/qiime_orchard_sequences.qza \
  --i-reference-taxonomy final_outputs/qiime_orchard_taxonomy.qza \
  --p-perc-identity 0.97 \
  --p-weak-id 0.94 \
  --p-query-cov 0.85 \
  --p-maxaccepts 12 \
  --p-maxrejects 256 \
  --p-min-consensus 0.8 \
  --p-threads 22 \
  --p-top-hits-only \
  --output-dir "$OUTDIR"
# results will be saved as classification.qza and search_results.qza in that dir


##############################################################################
#*               CHECK RESULTS
##############################################################################
#this code is kinda a mess... some of it is useful, some of it is not

#define out dir from above
OUTDIR="orchards/vsearch_out_20250912_205351"

# export taxonomy + hits
qiime tools export --input-path "$OUTDIR/classification.qza" --output-path "$OUTDIR/exported"
qiime tools export --input-path "$OUTDIR/search_results.qza" --output-path "$OUTDIR/exported"

# how deep did we classify? (species/genus/family/…/Unassigned)
awk -F'\t' 'NR>1{
  tax=$2
  if (tax=="Unassigned"){lvl="Unassigned"}
  else{
    n=split(tax,a,/;[ ]*/); lvl="kingdom"       # fallback
    for(i=n;i>=1;i--){
      if(a[i]!~/__NA$/){
        if(a[i]~/^s__/) lvl="species"
        else if(a[i]~/^g__/) lvl="genus"
        else if(a[i]~/^f__/) lvl="family"
        else if(a[i]~/^o__/) lvl="order"
        else if(a[i]~/^c__/) lvl="class"
        else if(a[i]~/^p__/) lvl="phylum"
        else if(a[i]~/^k__/) lvl="kingdom"
        break
      }
    }
  }
  c[lvl]++
} END{for(k in c) printf "%s\t%d\n", k, c[k] | "sort"}' "$OUTDIR/exported/taxonomy.tsv"


# top-hit %ID distribution 
awk 'NR==FNR{seen[$1]++; next} !q[$1]++ {print $1,$3}' \
  "$OUTDIR/exported/taxonomy.tsv" "$OUTDIR/exported/blast6.tsv" \
| awk '{id=int($2+0.5); c[id]++} END{for(i=80;i<=100;i++) if(c[i]) printf "%d%%\t%d\n", i, c[i]}' \
| sort -n

# get query lengths
qiime tools export --input-path orchards/rep-seqs_plate1.qza --output-path "$OUTDIR/rep_export"
awk '/^>/{id=substr($0,2);next}{len[id]+=length($0)}END{for(id in len)print id"\t"len[id]}' \
  "$OUTDIR/rep_export/dna-sequences.fasta" > "$OUTDIR/rep_export/qlen.tsv"

# top hit per query with %id and computed qcov
# blast6 cols: 1=qid 3=pident 7=qstart 8=qend
awk 'NR==FNR{L[$1]=$2;next} !seen[$1]++ {cov=($8-$7+1)/L[$1]; printf "%s\t%.2f\t%.3f\n",$1,$3,cov}' \
  "$OUTDIR/rep_export/qlen.tsv" "$OUTDIR/exported/blast6.tsv" > "$OUTDIR/exported/top_hits_id_cov.tsv"

# qcov distribution
awk '{cov=int($3*100+0.5); c[cov]++} END{for(i=50;i<=100;i++) if(c[i]) printf "%d%%\t%d\n", i, c[i]}' \
  "$OUTDIR/exported/top_hits_id_cov.tsv" | sort -n
# 100%    13865

##why were they unassigned?
# join top-hit metrics with assignment depth (0=Unassigned)
awk -F'\t' 'NR>1{r=($2=="Unassigned")?0:gsub(/;[ ]*/,";",$2)+1; print $1"\t"r}' \
  "$OUTDIR/exported/taxonomy.tsv" > "$OUTDIR/exported/assigned_depth.tsv"

join -t $'\t' -1 1 -2 1 \
  <(sort -k1,1 "$OUTDIR/exported/top_hits_id_cov.tsv") \
  <(sort -k1,1 "$OUTDIR/exported/assigned_depth.tsv") \
| awk -F'\t' '{id=$2; cov=int($3*100+0.5); d=$4; b=int(id/1); cc[b FS d]++}
              END{for(k in cc) print k, cc[k]}' \
| sort -k1,1n -k2,2n | column -t | sed '1i%ID_bin\tDepth(ranks)\tN'


# list of unassigned feature ids
awk -F'\t' 'NR>1 && $2=="Unassigned"{print $1}' "$OUTDIR/exported/taxonomy.tsv" \
  | sort > "$OUTDIR/exported/unassigned_ids.txt"

# which unassigned had at least one hit reported?
cut -f1 "$OUTDIR/exported/blast6.tsv" | sort -u > "$OUTDIR/exported/hit_ids.txt"

# counts
echo -n "unassigned total: "; wc -l < "$OUTDIR/exported/unassigned_ids.txt"
echo -n "unassigned WITH hits: "; comm -12 "$OUTDIR/exported/unassigned_ids.txt" "$OUTDIR/exported/hit_ids.txt" | wc -l
echo -n "unassigned with NO hits: "; comm -23 "$OUTDIR/exported/unassigned_ids.txt" "$OUTDIR/exported/hit_ids.txt" | wc -l


# VIZ
CLASS="$OUTDIR/classification.qza"
SEARCH="$OUTDIR/search_results.qza"

qiime metadata tabulate --m-input-file "$CLASS"  --o-visualization "$OUTDIR/classification.qzv" &
qiime metadata tabulate --m-input-file "$SEARCH" --o-visualization "$OUTDIR/search_results.qzv" &


##############################################################################
##############################################################################

######!    END OF CODE THAT I'VE RUN   ######

##############################################################################
##############################################################################




##############################################################################
#!        train naive bayes classifier (NOT USED, yet?)
#!  TOO MANY SEQs --> SWITCH TOVSEARCH CONSENSUS CLASSIFIER (SEE ABOVE)
##############################################################################
# memory logging
/usr/bin/time -v qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads final_outputs/qiime_orchard_sequences.qza \
  --i-reference-taxonomy final_outputs/qiime_orchard_taxonomy.qza \
  --o-classifier final_outputs/qiime_orchard_classifier.qza



##############################################################################
#*                3.  VSEARCH CONSENSUS CLASSIFIER 2 PASS METHOD
##############################################################################
#! DID NOT RUN YET:
# PASS 1: EXACT ONLY or very strict:
# ... --p-perc-identity 0.99 --p-query-cov 0.9 --p-maxaccepts 20 --o-classification pass1.qza
EXACT_OUT="orchards/vsearch_exact_$(date +%Y%m%d_%H%M%S)"

qiime feature-classifier classify-consensus-vsearch \
  --i-query orchards/rep-seqs_plate1.qza \
  --i-reference-reads final_outputs/qiime_orchard_sequences.qza \
  --i-reference-taxonomy final_outputs/qiime_orchard_taxonomy.qza \
  --p-search-exact \
  --p-threads 22 \
  --output-dir "$EXACT_OUT"

# PASS 2: more lenient to re-classify only the unassigned/high-rank reads with relaxed thresholds for genus/family:
# extract unassigned/high-rank reads, then:
#... --p-perc-identity 0.97 --p-query-cov 0.8 --p-maxaccepts 50 --o-classification pass2.qza
# merge the results;

# next:
    #--p-maxaccepts 50 (from 10) -->  better consensus voting at higher ranks.
    # --p-perc-identity 0.96 (or 0.95) --> may convert some “no consensus” cases to genus/family; keep an eye on species calls.
    # --p-min-consensus 0.7 more inclusive genus-level calls without lowering identity.



