# katherine carbeck
# primer trimming (single pass, anchored at 5' end of each read)
# using gene-specific ANML primers (without TruSeq tails)


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
# the goal is to determine how much we should truncate the reads before the paired end reads are joined. This will depend on the length of our amplicon, and the quality of the reads.
qiime demux summarize \
  --i-data trim_plate1.qza \
  --o-visualization trim_plate1.qzv

# How much of the total sequence do we need to preserve and still have a sufficient overlap to merge the paired end reads?
# How much of the poor quality sequence can we truncate before trying to merge?

