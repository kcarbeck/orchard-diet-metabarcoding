# Songbird diet metabarcoding: COI classifier (ANML) and diet metabarcoding pipeline

Build and apply a QIIME 2–compatible COI classifier for songbird diet metabarcoding using ANML primers. This repository provides modular shell scripts for: downloading COI references, importing and merging with CRABS, in silico PCR for ANML, dereplication and filtering, subsetting to project taxa, exporting QIIME 2 formats, training a Naive Bayes (QIIME 2 feature-classifier) model, and evaluating outputs.

The diet analysis workflow relies on QIIME 2 (2024.10) and R for `decontam` and posthoc analyses. Scripts live in `2_diet_analysis/` and are intended to be run in order, checking `.qzv` visualizations along the way.

Please site this repo if you use it! :) 

## Repository structure

```text
orchard-diet-metabarcoding/
├── 1_build_classifier/
│   ├── 00_setup.sh                    # create CRABS env,
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
│   ├── 09_evaluate_classifier.sh      # evaluate NB classifier on held-out set / metrics
│   ├── 10_confusion_matrix.sh         # build confusion matrix visualizations
│   └── helper_scripts/
│       ├── 06_helper_count_pests.sh   # presence/absence of priority pests in DB
│       └── gbif_filter_states.R       # GBIF geographic filters
├── 2_diet_analysis/
│   ├── 01_cutadapt.sh                 # primer trimming
│   ├── 02_dada2.sh                    # denoise paired reads
│   ├── 03_classify.sh                 # classify reads using trained classifier
│   ├── 04.1_decontam_prevalence.R     # identify contaminants (R)
│   ├── 04.2_filter_contaminants.sh    # remove contaminants per plate
│   ├── 05_merge_plates_and_years.sh   # merge taxonomy/metadata; merge artifacts
│   ├── 06_qc_and_get_analysis_files.sh# QC plots; analysis-only tables
│   ├── 07_posthoc_analyses.R          # alpha/beta; FOO/RRA; pest summaries (R)
│   ├── merge_taxonomy.R               # helper for taxonomy merge (R)
│   └── merge_metadata.R               # helper for metadata merge (R)
└── README.md
```

## Dependencies & install

Required software:

- QIIME 2: `qiime2-amplicon-2024.10` (activates via `source /programs/miniconda3/bin/activate qiime2-amplicon-2024.10` on Cornell BioHPC)
- CRABS: `>=1.8,<2.0` (banner shows v1.9.0 during help)
- cutadapt: `4.4` (required by CRABS in-silico PCR path)
- xopen: `<2.0` (workaround for known build issue)
- VSEARCH: `2.16.0` (used by CRABS pairwise alignment)
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



## Diet analysis pipeline (QIIME 2 + R)

### Expected inputs
- Per-plate QIIME 2 artifacts: `table_<plate>.qza`, `rep-seqs_<plate>.qza`, and `taxonomy_<plate>.qza` (after classification).
- Per-plate metadata TSVs with at least: `SampleID` (or `#SampleID`), `Species` (used for control detection), `Year`, `Site`, `Plate`.
- Negative controls labeled in `Species` as one of: `EBLANK`, `PBLANK`, `BLANK`, `EMPTY`; positive controls labeled `POS`.

### Run order
1. 01_cutadapt.sh
   - Primer trimming on demultiplexed reads; review `trim_*.qzv`.
2. 02_dada2.sh
   - Denoise paired-end reads; review `denoise_*.qzv` and `table_*.qzv`.
3. 03_classify.sh
   - Classify representative sequences with your trained classifier; review barplots.
4. 04.1_decontam_prevalence.R (Rscript)
   - Edit file paths at the top, then run.
   - Outputs: `decontam_*` directory with threshold sweep, figures, and `contaminant_feature_ids_thr_*.txt`.
5. 04.2_filter_contaminants.sh
   - Removes flagged ASVs from each per-plate table and filters rep-seqs to match.
6. 05_merge_plates_and_years.sh
   - Exports per-plate taxonomy; run `2_diet_analysis/merge_taxonomy.R` to create `taxonomy_merged.tsv`, then imports as `taxonomy_merged.qza`.
   - Runs `2_diet_analysis/merge_metadata.R` to produce `all_plates_metadata.tsv`; if duplicates are detected across plates, writes `rename_maps/*.tsv` and (optionally) renames per-plate sample IDs before merging.
   - Merges cleaned per-plate tables and rep-seqs into `table_merged_nocontam.qza` and `rep-seqs_merged_nocontam.qza`.
   - Validates metadata and checks coverage (reports any `missing_tax_ids.txt`).
7. 06_qc_and_get_analysis_files.sh
   - Generates QC barplots with controls.
   - Creates `table_analysis_only.qza` (samples only) and `table_analysis_only_arthropoda.qza` for downstream analyses.
   - Optional alpha-rarefaction visualization to guide depth selection (methods note only; do not use rarefied counts for DEICODE/ANCOM-BC/ANCOM-BC2).
8. 07_posthoc_analyses.R (Rscript)
   - Imports QZA artifacts directly into `phyloseq`.
   - Factors `Year`, `Plate`, and `Site` (orchard) if present.
   - Alpha metrics (Observed, Shannon) saved to `analysis/alpha_metrics.csv`.
   - Beta diversity: Jaccard and Bray–Curtis; PERMANOVA tests controlling for `Plate` and `Site`.
   - FOO and RRA summaries at Family/Genus; Year×Site CSVs and plots in `analysis/`.
   - Optional pest layer: set `pest_csv` (default `data/pest_taxa.csv` with columns `rank` ∈ {Genus, Species}, `taxon`). Outputs per-sample and Year×Site pest summaries and top pest taxa tables in `analysis/`.
   - Example differential abundance at Family level via ANCOM-BC2 (saved in `analysis/`).

### Notes and tips
- Keep `taxonomy_merged.qza`, `rep-seqs_merged_nocontam.qza`, and `all_plates_metadata.tsv` as your canonical merged artifacts for provenance.
- Use unrarefied counts for DEICODE/RPCA, ANCOM‑BC/ANCOM‑BC2, and FOO/RRA summaries. Use rarefaction only for visualization or where even depth is explicitly required.
- The `Site` metadata column is used as orchard ID in posthoc analyses and figures.

### Optional extras
- `misc_scripts/07_diversity.sh` provides QIIME-native alpha/beta pipelines (including DEICODE) and exports grouped FOO/RRA tables for R plotting.

## Build classifier pipeline (Naive Bayes)

Scripts in `1_build_classifier/` train and evaluate a QIIME 2 Naive Bayes classifier tailored to the ANML COI amplicon.

1. 00_setup.sh
   - Create CRABS environment and ensure QIIME 2 is available.
2. 01_download_COI.sh
   - Download raw COI references from BOLD/NCBI/MIDORI via CRABS.
3. 02_merge.sh
   - Import, merge, and deduplicate references.
4. 03_in_silico_pcr.sh
   - In-silico PCR with ANML primers to isolate the target amplicon region.
5. 04_global_alignment.sh
   - Global alignment (VSEARCH) for length QC and amplicon validation.
6. 05_database_filtering.sh
   - Dereplicate to unique species; quality and length filters.
7. 06.1_database_subsetting.sh / 06.2_gbif_subsetting.sh
   - Optional biological/geographic subsetting to project-relevant taxa.
8. 07.1_export.sh and 07.2_clean_db.sh
   - Export QIIME 2 `FeatureData[Sequence]` and `FeatureData[Taxonomy]`; optional taxonomy cleanup.
9. 08_train_nb_classifier.sh
   - Train Naive Bayes classifier (QIIME 2 `feature-classifier classify-sklearn` compatible); ensure primer-trimmed amplicon sequences and matching taxonomy are used.
10. 09_evaluate_classifier.sh
   - Evaluate classifier against a held-out or benchmark set; export accuracy/precision/recall and barplots.
11. 10_confusion_matrix.sh
   - Generate confusion matrices across taxonomic ranks for visual QA.

Notes:
- Ensure the training sequences match the exact ANML amplicon region used for classification (primer-trimmed, consistent orientation).
- Prefer species-level dereplication for cleaner NB training; retain lineage strings with consistent rank prefixes.




## Citations

- QIIME 2: Bolyen E, Rideout JR, et al. (2019) Nature Biotechnology 37, 852–857. [Project page](https://qiime2.org)
- VSEARCH: Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) PeerJ 4:e2584. [Project page](https://github.com/torognes/vsearch)
- CRABS: Jeunen G-J, Dowle E, Edgecombe J, von Ammon U, Gemmell NJ, Cross H. (2022) Molecular Ecology Resources. doi:10.1111/1755-0998.13741. [Docs](https://github.com/GenomicsAotearoa/crabs)
- ANML primers: [REFERENCE FOR PRIMERS HERE]

