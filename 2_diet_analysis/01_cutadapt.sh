# primer trimming (single pass, anchored at 5' end of each read)
# using gene-specific ANML primers (without TruSeq tails)

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source /programs/miniconda3/bin/activate qiime2-amplicon-2024.10

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux_plate1.qza \
  --p-front-f    GGTCAACAAATCATAAAGATATTGG \
  --p-front-r    GGWACTAATCAATTTCCAAATCC \
  --p-adapter-f  GGATTTGGAAATTGATTAGTWCC \
  --p-adapter-r  CCAATATCTTTATGATTTGTTGACC \
  --p-match-adapter-wildcards \
  --p-match-read-wildcards \
  --p-cores 22 \
  --o-trimmed-sequences trimmed_plate1.qza \
  --verbose > cutadapt_out_plate1.txt


# inspect .qzv after trimming to see where quality tapers off
qiime demux summarize \
  --i-data trim_plate1.qza \
  --o-visualization trim_plate1.qzv
