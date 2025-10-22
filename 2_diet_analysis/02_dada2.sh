

# suggested params:
# - expected overlap ≈ truncF + truncR − insert ≈ 130 + 130 − 181 = 79 bp (plenty).
# - --p-min-overlap 20 gives headroom if you ever need to trim a bit shorter for quality.

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trim_plate1.qza \
  --p-trunc-len-f 125 \
  --p-trunc-len-r 125 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-min-overlap 20 \
  --p-n-threads 22 \
  --o-representative-sequences rep-seqs_plate1.qza \
  --o-table table_plate1.qza \
  --o-denoising-stats denoise_plate1.qza

qiime metadata tabulate \
  --m-input-file denoise_plate1.qza \
  --o-visualization denoise_plate1.qzv

qiime feature-table summarize \
  --i-table table_plate1.qza \
  --m-sample-metadata-file 2024_wrangled_data_2_CLEAN.tsv \
  --o-visualization table_plate1.qzv
