# katherine carbeck
# import and demultiplex reads

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source /programs/miniconda3/bin/activate qiime2-amplicon-2024.10

qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /path/to/your/data  \
  --output-path demux_plate1.qza

qiime demux summarize \
 --i-data demux_plate1.qza \
 --o-visualization demux_plate1.qzv

