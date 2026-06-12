# katherine carbeck
# dada2 denoising


#DADA2 denoising models run-specific, quality-score–dependent sequencing errors to statistically infer the true biological sequences in your samples (ASVs). It de-replicates reads, learns an error model, and tests whether each unique read could be an error of a more abundant one—if not, it’s retained as its own ASV. During denoising it also merges paired-end reads (using overlap) and removes chimeras. because error profiles differ by run, denoise each sequencing run separately in qiime2, then merge the resulting feature tables/rep-seqs afterward.


# suggested params:
# - expected overlap ≈ truncF + truncR − insert ≈ 130 + 130 − 181 = 79 bp (plenty).
# - --p-min-overlap 20 gives headroom if you ever need to trim a bit shorter for quality.

# our ANML amplicon is ~180 bp between our primers. Our goal overlap is ~50bp. trunc 125/125 is good --> overlap=250-180=70bp which is plenty. If forward-favored and R2 is weaker can adjust to something like 140/110 (250-180=70bp). 
# Could bump min overlap up if we want to be more conservative.
# many conservative ANML pipelines use 130/130 as trunc lengths, but 135/125 looked good for us. 


##############################################################################
#*                FILE PATHS (UPDATE THESE ACCORDINGLY)
##############################################################################
WORK_DIR=/lustre2/home/lc736_0001/orchards/
DIR_2024=$WORK_DIR/orchards_diet_2024
DIR_2025=$WORK_DIR/orchards_diet_2025

# 2024
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $DIR_2024/trimmed_plate1_2024.qza \
  --p-trunc-len-f 125 \
  --p-trunc-len-r 125 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-min-overlap 20 \
  --p-n-threads 22 \
  --o-representative-sequences $DIR_2024/rep-seqs_plate1_2024.qza \
  --o-table $DIR_2024/table_plate1_2024.qza \
  --o-denoising-stats $DIR_2024/denoise_plate1_2024.qza &














qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trim_plate1.qza \
  --p-trunc-len-f 130 \
  --p-trunc-len-r 130 \
  --p-trim-left-f 1 \
  --p-trim-left-r 1 \
  --p-min-overlap 20 \
  --p-n-threads 22 \
  --o-representative-sequences rep-seqs_plate1.qza \
  --o-table table_plate1.qza \
  --o-denoising-stats denoise_plate1.qza

# metadata on denoising
qiime metadata tabulate \
  --m-input-file denoise_plate1.qza \
  --o-visualization denoise_plate1.qzv
# look at denoise_plate1.qzv to see how many seqs passed the cutoff for each sample at each step of the denoising process. 

# table of per-sample sequence counts
qiime feature-table summarize \
  --i-table table_plate1.qza \
  --m-sample-metadata-file 2024_wrangled_data_2_CLEAN.tsv \
  --o-visualization table_plate1.qzv

# unique sequences accross samples
qiime feature-table tabulate-seqs\
   --i-data rep-seqs_plate1.qza\
   --o-visualization rep-seqs_plate1.qzv
# see sequences and the distribution of seq lengths. each seq is a link to blast against ncbi


##############################################################################
#* export sequences to fasta and check distribution of sequence lengths
##############################################################################

qiime tools export \
  --input-path /lustre2/home/lc736_0001/orchards/2025/2025/rep-seqs_plate2.qza \
  --output-path /lustre2/home/lc736_0001/orchards/orchards_diet_2025/repseqs_export_2025_plate2

# cd to the export directory and check the distribution of sequence lengths
awk '/^>/{next}{print length($0)}' dna-sequences.fasta \
| sort -n | awk '{a[++n]=$1} END{printf("min=%d p1=%d median=%d p99=%d max=%d\n",a[1],a[int(0.01*n)],a[int(0.5*n)],a[int(0.99*n)],a[n])}'

awk '/^>/{next}{print length($0)}' dna-sequences.fasta \
  | sort -n | uniq -c
