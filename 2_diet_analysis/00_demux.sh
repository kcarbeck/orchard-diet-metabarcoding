# katherine carbeck
# import and demultiplex reads

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source /programs/miniconda3/bin/activate qiime2-amplicon-2024.10


qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_plate1_2024.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path /lustre2/home/lc736_0001/orchards/2024/demux_plate1.qza

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_plate1_2025.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path /lustre2/home/lc736_0001/orchards/2025/demux_plate1.qza


qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_plate2_2025.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path /lustre2/home/lc736_0001/orchards/2025/demux_plate2.qza


##############################################################################
#*                 summarize demultiplexed reads
##############################################################################

qiime demux summarize \
 --i-data demux_plate1.qza \
 --o-visualization demux_plate1.qzv &

qiime demux summarize \
 --i-data demux_plate1.qza \
 --o-visualization demux_plate1.qzv &

qiime demux summarize \
 --i-data demux_plate2.qza \
 --o-visualization demux_plate2.qzv &

