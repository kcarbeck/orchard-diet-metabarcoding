# Songbird diet metabarcoding: COI classifier (ANML) build-and-use

Build and apply a QIIME 2–compatible COI classifier for songbird diet metabarcoding using the ANML primers. This repository provides modular shell scripts (no workflow manager) for: downloading COI references, importing and merging with CRABS, in silico PCR for ANML, dereplication and filtering, subsetting to project taxa, exporting QIIME 2 formats, training or using a classifier (VSEARCH consensus recommended at this scale), and evaluating outputs.

Key tools: CRABS (reference curation), cutadapt (via CRABS), VSEARCH, QIIME 2 (2024.10). Utilities include GNU parallel, awk, and basic R for comparisons.

## Repository structure

```text
orchard-diet-metabarcoding/
├── 1_build_classifier/
│   ├── 00_setup.sh                    # create CRABS env,
│   ├── 01_download_COI.sh             # CRABS downloads from BOLD, NCBI, MIDORI
│   ├── 02_merge.sh                    # CRABS import; merge + dedup
│   ├── 03_in_silico_pcr.sh            # clean inputs; ANML in-silico PCR 
│   ├── 04_global_alignment.sh         # recover amplicons (VSEARCH pairwise); length figure
│   ├── 05_database_filtering.sh       # dereplicate (unique_species) + quality filter
│   ├── 06_database_subsetting.sh      # subset database to relevant taxa
│   ├── 07_export.sh                   # export QIIME2 sequences/taxonomy (SINTAX and BLAST-tax if needed)
│   ├── 08_train_qiime_classifier.sh   # VSEARCH classification
│   ├── 09_evaluate.sh                 # export classification/search; QC metrics; barplots
│   └── helper_scripts/
│       ├── 06_helper_count_pests.sh   # presence/absence of priority pests in DB
│       └── 06A_TEST_SUBSET.sh         # tiny test subset pipeline from in silico PCR output
├── 2_diet_analysis/
│   ├── 01_import_demux.sh             # importing demultiplexed reads 
│   ├── 02_cutadapt.sh                 # primer trimming
│   └── 03_dada2.sh                    # DADA2 
│   └── ......                         # fill in pipeline
└── README.md
```

## Dependencies & install

Required software:

- QIIME 2: `qiime2-amplicon-2024.10` (activates via `source /programs/miniconda3/bin/activate qiime2-amplicon-2024.10` on Cornell BioHPC)
- CRABS: `>=1.8,<2.0` (banner shows v1.9.0 during help)
- cutadapt: `4.4` (required by CRABS in-silico PCR path)
- xopen: `<2.0` (workaround for known build issue)
- VSEARCH: `2.16.0` (used by CRABS pairwise alignment and QIIME’s vsearch)
- BLAST+ makeblastdb: `2.10.1+` (listed as external dependency)
- GNU parallel, awk; optional R (base) for comparisons

Create and populate CRABS environment (Conda/Mamba):

```bash
conda create -y -n CRABS python=3.10
conda activate CRABS
conda config --env --set solver classic
conda config --env --set channel_priority strict
conda install -y -c conda-forge -c bioconda \
  python=3.10 \
  "xopen<2.0" \
  "cutadapt>=4.3,<5.0" \
  "biopython>=1.80" \
  "crabs>=1.8,<2.0"

# verify
export PATH="$CONDA_PREFIX/bin:$PATH" && hash -r
type -a crabs
crabs --help | head
```

QIIME 2 activation (Cornell BioHPC example):

```bash
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
source /programs/miniconda3/bin/activate qiime2-amplicon-2024.10
qiime --help | head
```


## Citations

- QIIME 2: Bolyen E, Rideout JR, et al. (2019) Nature Biotechnology 37, 852–857. [Project page](https://qiime2.org)
- VSEARCH: Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) PeerJ 4:e2584. [Project page](https://github.com/torognes/vsearch)
- CRABS: Jeunen G-J, Dowle E, Edgecombe J, von Ammon U, Gemmell NJ, Cross H. (2022) Molecular Ecology Resources. doi:10.1111/1755-0998.13741. [Docs](https://github.com/GenomicsAotearoa/crabs)
- ANML primers: [REFERENCE FOR PRIMERS HERE]

