# setup env and dirs for metabarcoding pipeline
# katherine carbeck
# 07 sep 2025

# CITE: 
# Jeunen, G.-J., Dowle, E., Edgecombe, J., von Ammon, U., Gemmell, N. J., & Cross, H. (2022). crabs—A software program to generate curated reference databases for metabarcoding sequencing data. Molecular Ecology Resources, 00, 1– 14.](https://doi.org/10.1111/1755-0998.13741)

# required dependencies:
#python modules:
   # requests v 2.32.3 - A Python library for making HTTP requests.
   # rich v 13.3.5 - A Python library for beautiful terminal formatting.
   # rich-click v 1.7.2 - A Python library that extends Click to provide rich formatting for command-line interfaces.
   # matplotlib v 3.8.0 - A plotting library used for creating static, animated, and interactive visualizations.
   # numpy v 1.26.2 - A library for working with arrays and numerical operations.

#external software programs:
   # makeblastdb v 2.10.1+ - A tool from the BLAST+ suite for creating BLAST databases.
   # cutadapt v 4.4 - A tool for trimming adapter sequences from high-throughput sequencing reads.
   # vsearch v 2.16.0 - A tool for metagenomics data analysis and high-throughput sequencing data analysis.
   # clustalw2 v 2.1 - A tool for multiple sequence alignments.
   # FastTree v 2.1.11 - A tool for inferring approximately maximum-likelihood phylogenetic trees.



conda activate CRABS
# use the classic solver just for this env (avoids the libmamba warning storm)
conda config --env --set solver classic
conda config --env --set channel_priority strict

# install with pins that avoid the bad xopen build error
conda install -c conda-forge -c bioconda \
  python=3.10 \
  "xopen<2.0" \
  "cutadapt>=4.3,<5.0" \
  "biopython>=1.80" \
  "crabs>=1.8,<2.0"

# check if worked
crabs --version
which crabs
crabs --help | head


which crabs
/lustre2/home/lc736_0001/diet/songbird_coi_database/reference_database_creator/crabs
# that path is shadowing the conda one. also, crabs doesnt support --version (the banner already shows v1.9.0), and those --exclude warnings are harmless.


# make sure the conda env's bin comes first
export PATH="$CONDA_PREFIX/bin:$PATH"
hash -r   # clear any hashed path to 'crabs'

# confirm we’re hitting the conda one
type -a crabs
# should show: /home/kcarbeck/miniconda3/envs/CRABS/bin/crabs ... at the top

# sanity-check
crabs --help | head





##############################################################################
#*                 the manual download version
##############################################################################
#export path
export PATH="/lustre2/home/lc736_0001/diet/songbird_coi_database/reference_database_creator:$PATH"

# direct path to crabs
crabs --help
