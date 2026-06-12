# Songbird diet metabarcoding: COI classifier (ANML) and diet metabarcoding pipeline

Build and apply a QIIME 2–compatible COI classifier for songbird diet metabarcoding using ANML primers. This repository provides modular shell scripts for: downloading COI references, importing and merging with CRABS, in silico PCR for ANML, dereplication and filtering, subsetting to project taxa, exporting QIIME 2 formats, training a Naive Bayes (QIIME 2 feature-classifier) model, and evaluating outputs.

The diet analysis workflow relies on QIIME 2 (2024.10) and R for `decontam` and posthoc analyses. Scripts live in `2_diet_analysis/` and are intended to be run in order, checking `.qzv` visualizations along the way.

**Please cite this repo if you use it! :)**

## Repository structure

```text
orchard-diet-metabarcoding/
├── 1_build_classifier/
│   ├── 00_setup.sh                    # create CRABS env
│   ├── 01_download_COI.sh             # CRABS downloads from BOLD, NCBI, MIDORI
│   ├── 02_merge.sh                    # CRABS import; merge + dedup
│   ├── 03_in_silico_pcr.sh            # ANML in-silico PCR
│   ├── 04_global_alignment.sh         # recover amplicons; length QC
│   ├── 05_database_filtering.sh       # dereplicate + filter
│   ├── 06.1_database_subsetting.sh    # subset database to target taxa
│   ├── 06.2_gbif_subsetting.sh        # GBIF-based subsetting (R helper)
│   ├── 07.1_export.sh                 # export QIIME 2 sequences/taxonomy
│   ├── 07.2_clean_db.sh               # taxonomy cleanup
│   ├── 08_train_nb_classifier.sh      # train Naive Bayes classifier (QIIME 2)
│   ├── 09_evaluate_classifier.sh      # evaluate NB classifier on held-out set
│   ├── 10_confusion_matrix.sh         # build confusion matrix visualizations
│   └── helper_scripts/
│       ├── 06_helper_count_pests.sh   # presence/absence of priority pests in DB
│       └── gbif_filter_states.R       # GBIF geographic filters
├── 2_diet_analysis/
│   ├── 00_demux.sh                    # import and demultiplex raw reads
│   ├── 01_cutadapt.sh                 # primer trimming
│   ├── 02_dada2.sh                    # denoise paired reads
│   ├── 03_classify.sh                 # classify reads using trained classifier
│   ├── 04.1_merge_metadata.R          # merge metadata across plates/years
│   ├── 04.2_merge_plates_and_years.sh # merge tables, rep-seqs, taxonomy
│   ├── 05_decontam_prevalence.R       # identify contaminants using negatives
│   ├── 06_qc_and_get_analysis_files.sh# quality filtering, rarefaction curves, SRS
│   ├── 07_posthoc_analyses.R          # alpha/beta; FOO/RRA; pest summaries
│   ├── helper_scripts/
│   │   ├── merge_taxonomy.R           # helper for taxonomy merge
│   │   ├── 03.1_evaluate_classifier.sh
│   │   └── 03.2_confusion_matrix.sh
│   └── README.md                      # detailed diet pipeline documentation
└── README.md
```

## Dependencies & install

Required software:

- QIIME 2: `qiime2-amplicon-2024.10` (activates via `source /programs/miniconda3/bin/activate qiime2-amplicon-2024.10` on Cornell BioHPC)
- CRABS: `>=1.8,<2.0` 
- cutadapt: `4.4`
- xopen: `<2.0` 
- VSEARCH: `2.16.0` 
- BLAST+ makeblastdb: `2.10.1+` 
- GNU parallel, awk; R (base) 

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

## Diet analysis pipeline (QIIME 2 + R)

### Pipeline overview

```
Raw Reads → Demux → Cutadapt → DADA2 → Classify → Merge → Decontam → QC/Filter → Analysis
```

### Expected inputs
- Per-plate demultiplexed reads (from sequencing facility)
- Per-plate metadata TSVs with: `SampleID` (or `#SampleID`), `Species` (for control detection), `Year`, `Site`, `PlateName`
- Negative controls labeled in `Species` as: `EBLANK`, `PBLANK`, `BLANK`, `EMPTY`
- Positive controls labeled as: `POS`

### Run order

#### Steps 00-03: Run for each plate separately

1. **00_demux.sh** - Import and demultiplex raw sequencing reads
2. **01_cutadapt.sh** - Primer trimming; review `trim_*.qzv`
3. **02_dada2.sh** - Denoise paired-end reads; review `denoise_*.qzv` and `table_*.qzv`
4. **03_classify.sh** - Classify representative sequences with trained classifier; review barplots

#### Steps 04+: Run once after all plates are processed

5. **04.1_merge_metadata.R** (Rscript)
   - Merges metadata from all plates/years
   - Creates `Year_Plate` column for batch identification in decontam
   - Prefixes duplicate control IDs to make them unique
   - Outputs: `all_plates_metadata.tsv`, `rename_maps/*.tsv` (if duplicates found)

6. **04.2_merge_plates_and_years.sh**
   - Merges feature tables, rep-seqs, and taxonomy from all plates
   - Calls `helper_scripts/merge_taxonomy.R` to merge taxonomy files
   - Outputs: `merged_output/table_merged.qza`, `rep-seqs_merged.qza`, `taxonomy_merged.qza`

7. **05_decontam_prevalence.R** (Rscript)
   - Edit config block with file paths, then run
   - Uses prevalence method with batch mode (Year_Plate column)
   - Outputs: `decontam_output/` with threshold sweep, diagnostic plots, `contaminant_feature_ids_thr_*.txt`
   - After running, filter contaminants using commands at end of 04.2 script

8. **06_qc_and_get_analysis_files.sh**
   - Quality filtering (removes low-depth samples and rare ASVs)
   - Generates rarefaction curves (diagnostic only - does NOT rarefy data)
   - Applies SRS normalization for alpha diversity
   - Outputs: `analysis_tables/table_arthropoda.qza` (unrarefied), `table_arthropoda_srs.qza` (SRS-normalized)

9. **07_posthoc_analyses.R** (Rscript)
   - Imports QZA artifacts directly into `phyloseq`
   - Alpha diversity metrics (Observed, Shannon, Chao1)
   - Beta diversity (Jaccard, Bray-Curtis) with PERMANOVA
   - FOO and RRA summaries at Family level
   - Optional pest analysis and ANCOM-BC2 differential abundance
   - Outputs: `analysis/` directory with CSVs and plots

### Key terminology: rarefaction vs rarefying

| Term | What it is | Do we use it? |
|------|------------|---------------|
| **Rarefaction curves** | Diagnostic visualization to assess sampling depth | YES (step 06) |
| **Rarefying (subsampling)** | Data transformation to equalize depth | NO - use SRS instead |

### Which table to use for each analysis

| Analysis | Table | Rationale |
|----------|-------|-----------|
| Rarefaction curves | Unrarefied | Diagnostic only |
| Alpha diversity | SRS-normalized | Depth-sensitive metric |
| Beta diversity (PERMANOVA) | Unrarefied | Bray-Curtis handles uneven depth |
| FOO/RRA | Unrarefied | Presence/absence less depth-sensitive |
| Differential abundance | Unrarefied | ANCOM-BC has internal normalization |

### Quality filtering thresholds

Default thresholds in `06_qc_and_get_analysis_files.sh`:

| Filter | Default | Description |
|--------|---------|-------------|
| `MIN_SAMPLE_READS` | 1000 | Remove samples with fewer reads |
| `MIN_ASV_READS` | 10 | Remove ASVs with fewer total reads |
| `MIN_ASV_SAMPLES` | 2 | Remove ASVs appearing in only 1 sample |

Adjust based on your data and rarefaction curves.

### Notes and tips
- Keep `taxonomy_merged.qza`, `rep-seqs_merged.qza`, and `all_plates_metadata.tsv` as your canonical merged artifacts for provenance
- The `Year_Plate` column is used for batch mode in decontam (accounts for plate-to-plate variation)
- The `Site` metadata column is used as orchard ID in posthoc analyses


## Build classifier pipeline (Naive Bayes)

Scripts in `1_build_classifier/` train and evaluate a QIIME 2 Naive Bayes classifier tailored to the ANML COI amplicon.

1. 00_setup.sh - Create CRABS environment and ensure QIIME 2 is available
2. 01_download_COI.sh - Download raw COI references from BOLD/NCBI/MIDORI via CRABS
3. 02_merge.sh - Import, merge, and deduplicate references
4. 03_in_silico_pcr.sh - In-silico PCR with ANML primers to isolate target amplicon region
5. 04_global_alignment.sh - Global alignment (VSEARCH) for length QC and amplicon validation
6. 05_database_filtering.sh - Dereplicate to unique species; quality and length filters
7. 06.1_database_subsetting.sh / 06.2_gbif_subsetting.sh - Optional biological/geographic subsetting
8. 07.1_export.sh and 07.2_clean_db.sh - Export QIIME 2 formats; optional taxonomy cleanup
9. 08_train_nb_classifier.sh - Train Naive Bayes classifier
10. 09_evaluate_classifier.sh - Evaluate classifier against held-out or benchmark set
11. 10_confusion_matrix.sh - Generate confusion matrices across taxonomic ranks

Notes:
- Ensure training sequences match the exact ANML amplicon region (primer-trimmed, consistent orientation)
- Prefer species-level dereplication for cleaner NB training


## Citations

- QIIME 2: Bolyen E, Rideout JR, et al. (2019) Nature Biotechnology 37, 852–857. [Project page](https://qiime2.org)
- VSEARCH: Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) PeerJ 4:e2584. [Project page](https://github.com/torognes/vsearch)
- CRABS: Jeunen G-J, Dowle E, Edgecombe J, von Ammon U, Gemmell NJ, Cross H. (2022) Molecular Ecology Resources. doi:10.1111/1755-0998.13741. [Docs](https://github.com/GenomicsAotearoa/crabs)
- decontam: Davis NM, Proctor DM, et al. (2018) Microbiome 6:226. doi:10.1186/s40168-018-0605-2
- SRS: Beule L, Karlovsky P. (2020) PeerJ 8:e9593. doi:10.7717/peerj.9593
- Diet metabarcoding: Deagle BE, Thomas AC, et al. (2019) Molecular Ecology 28:1542-1558. doi:10.1111/mec.14734
