

# suggested params:
# - expected overlap ≈ truncF + truncR − insert ≈ 130 + 130 − 181 = 79 bp (plenty).
# - --p-min-overlap 20 gives headroom if you ever need to trim a bit shorter for quality.



qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed_plateX.qza \
  --p-trunc-len-f 130 \
  --p-trunc-len-r 130 \
  --p-min-overlap 20 \
  --p-trim-left-f 0 --p-trim-left-r 0 \
  --o-representative-sequences rep-seqs_plateX.qza \
  --o-table table_plateX.qza \
  --o-denoising-stats denoise_plateX.qza
