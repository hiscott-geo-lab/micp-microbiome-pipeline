# MICP Microbiome Analysis Pipeline

Reproducible R pipeline for microbial community analysis of MICP-treated soils. Developed for the ERDC-funded project *"Engineering Bio-Mediated Soil Improvement for Planted Infrastructure and Ecosystem Resiliency"* (NC State University, Montoya Geoenvironmental Lab).

Sequencing data processed via Nanopore MinION → NanoClust/NanoScript → T-BAS phylogenetic placement (LIFE1 reference tree). This pipeline picks up from T-BAS output and carries through all downstream community analyses.

---

## Quick Start

```r
# 1. Clone the repo and open RStudio, set working directory to repo root
setwd("/path/to/micp-microbiome-pipeline")

# 2. Place your input files in data/raw/ (see Input Files below)

# 3. Edit R/00_config.R for your experiment

# 4. Run the full pipeline
source("run_pipeline.R")
```

---

## Repository Structure

```
micp-microbiome-pipeline/
├── run_pipeline.R              # Master runner — start here
├── R/
│   ├── 00_config.R             # ← EDIT THIS for each experiment
│   ├── 01_build_phyloseq.R     # Data import, QC, phyloseq construction
│   ├── 02_relative_abundance.R # Stacked bar plots (genus-level)
│   ├── 03_alpha_diversity.R    # Shannon, Simpson, Observed, Faith's PD
│   ├── 04_beta_diversity.R     # PCoA, NMDS ordination
│   ├── 05_permanova.R          # PERMANOVA, pairwise, PERMDISP
│   └── 06_ncycling_guilds.R    # N-cycling functional guild analysis
├── data/
│   ├── raw/                    # Your input files go here (not tracked by git)
│   └── processed/              # ps.rds, ps_rarefied.rds (auto-generated)
├── outputs/
│   ├── figures/                # All plots (.tiff LZW + .png)
│   └── tables/                 # All statistical tables (.csv)
└── docs/
    └── methods_template.md     # Methods text template for manuscripts
```

---

## Input Files

Place these four files in `data/raw/` before running. File paths are set in `R/00_config.R`.

| File | Format | Description |
|---|---|---|
| `phyloseq_otu_table.csv` | CSV | T-BAS output: samples as rows, OTUs as columns, first column = `SampleID` |
| `phyloseq_taxonomy_table.csv` | CSV | T-BAS output: `FeatureID` + one column per rank (Kingdom through Genus) |
| `ERDC_Metadata_corrected.csv` | CSV | Sample metadata: `Sample_Code`, `Treatment`, `Species`, plus any covariates |
| `annotated_tree.nwk` | Newick | Phylogenetic tree from T-BAS placement. Optional — set `PATHS$tree <- NULL` in config to skip UniFrac and Faith's PD |

**OTU table note:** T-BAS exports samples as rows. The pipeline transposes this automatically to the taxa-as-rows convention required by phyloseq.

---

## Configuring for a New Experiment

All experiment-specific settings live in **`R/00_config.R`**. You should not need to edit any other file.

### Minimum changes for Task 2 / future datasets:

```r
# 1. Update experiment ID (used in filenames)
EXPERIMENT_ID <- "ERDC_Task2"

# 2. Update file paths
PATHS <- list(
  otu_table = "data/raw/task2_otu_table.csv",
  taxonomy  = "data/raw/task2_taxonomy_table.csv",
  metadata  = "data/raw/task2_metadata.csv",
  tree      = "data/raw/task2_tree.nwk"
)

# 3. Update treatment labels if recipe changed
TREATMENT_MAP <- list(
  raw_levels     = c("Control", "MICP 1"),
  display_labels = c("Control", "U250")
)

# 4. Update excluded samples for new QC decisions
EXCLUDE_SAMPLES <- c("BAD_SAMPLE_1")

# 5. Update PERMANOVA formula if metadata structure changes
# (e.g., Task 2 has plant richness and species composition as covariates)
PERMANOVA_FORMULA <- "~ Treatment + PlantRichness"
```

### Adding new plant species:
```r
# In ADAPTATION_GROUPS, add the new species code:
ADAPTATION_GROUPS <- c(
  ADAPTATION_GROUPS,    # keep existing entries
  NEWSP = "Dual tolerant"
)
```

### Adding or modifying N-cycling guilds:
```r
# In N_GUILDS, add genera to existing guilds or define a new guild:
N_GUILDS$Anammox <- c("Candidatus Kuenenia", "Candidatus Brocadia",
                       "Candidatus Jettenia")
```

---

## Outputs

All outputs are written to `outputs/` automatically.

### Figures (`outputs/figures/`)

| File | Description |
|---|---|
| `bar_treatment.tiff/.png` | Relative abundance by treatment, top N genera |
| `bar_by_species.tiff/.png` | Relative abundance faceted by species × treatment |
| `alpha_diversity.tiff/.png` | Shannon, Simpson, Observed OTUs by treatment |
| `alpha_faiths_pd.tiff/.png` | Faith's PD by treatment |
| `alpha_faiths_pd_by_species.tiff/.png` | Faith's PD faceted by species |
| `beta_braycurtis_pcoa.tiff/.png` | Bray-Curtis PCoA |
| `beta_braycurtis_nmds.tiff/.png` | Bray-Curtis NMDS |
| `beta_unifrac_unweighted_pcoa.tiff/.png` | Unweighted UniFrac PCoA |
| `beta_unifrac_weighted_pcoa.tiff/.png` | Weighted UniFrac PCoA |
| `betadisper_boxplot.tiff/.png` | Distance-to-centroid boxplot (PERMDISP) |
| `fig_ncycling_guilds_by_treatment.tiff/.png` | N-cycling guild relative abundance |

All TIFF figures use LZW compression at 300 DPI (journal submission standard).

### Tables (`outputs/tables/`)

| File | Description |
|---|---|
| `alpha_diversity_stats.csv` | Kruskal-Wallis results for alpha metrics |
| `alpha_faiths_pd_stats.csv` | Faith's PD summary by treatment |
| `permanova_results.csv` | PERMANOVA (all distance metrics) |
| `pairwise_permanova_bray.csv` | Pairwise PERMANOVA with Bonferroni correction |
| `pairwise_permanova_unifrac.csv` | " (Unweighted UniFrac) |
| `pairwise_permanova_wunifrac.csv` | " (Weighted UniFrac) |
| `ncycling_guild_summary.csv` | Guild-level relative abundance per sample |
| `ncycling_guild_stats_kw.csv` | Kruskal-Wallis results per guild |
| `ncycling_guild_stats_dunn.csv` | Pairwise Dunn tests per guild |

---

## Dependencies

Install once in R:

```r
install.packages(c(
  "tidyverse", "vegan", "patchwork", "ggrepel",
  "picante", "dunn.test", "rstatix", "ape"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
```

Tested with R ≥ 4.2. Package versions used in original analysis are listed in `docs/session_info.txt` (generated at end of a pipeline run).

---

## Statistical Approach Notes

### PERMANOVA
`adonis2(by = "terms")` is used throughout (vegan ≥ 2.6). Predictors are tested sequentially; order is set in `PERMANOVA_FORMULA` in config. The formula default `~ Treatment + Species` reflects the study design of Task 1 (one species per pot). For Task 2 (mixed species pots), update the formula to include `PlantRichness` or species composition covariates as appropriate.

### PERMDISP
Betadisper/PERMDISP is run as a validity check for PERMANOVA, not as a standalone ecological result. A significant PERMDISP indicates unequal within-group dispersions, which can inflate PERMANOVA significance. Report PERMDISP alongside PERMANOVA and interpret accordingly.

### N-cycling guilds
Guild membership is non-exclusive. Genera with documented roles in multiple N-transformation pathways (e.g., *Paenibacillus*) contribute their relative abundance independently to each applicable guild. Guild-level abundances do not sum to unity across guilds and should not be presented as if they do.

### Rarefaction
Beta diversity and PERMANOVA use a rarefied phyloseq object (`ps_rarefied.rds`). Alpha diversity is computed on the unrarefied object. The rarefaction depth defaults to the minimum sample depth; set `RAREFY_DEPTH` explicitly in config when one or two outlier-low samples would otherwise force a very shallow rarefaction that discards most data.

---

## Data Provenance

Raw sequencing data → NanoClust/NanoScript (Nextflow) → T-BAS v2.4 (LIFE1 reference tree) → this pipeline.

Bioinformatics: Dr. Ignazio Carbone (CIFR, NC State).  
Sequencing: Dr. Mary Anna Carbone (NC State).  
Pipeline development: Hannah Hiscott (NC State, Montoya Lab).

---

## Citation

If you use this pipeline, please cite the associated dissertation:

> Hiscott, H.F. (2026). *Engineering Bio-Mediated Soil Improvement for Planted Infrastructure and Ecosystem Resiliency.* PhD Dissertation, North Carolina State University.

And the relevant journal papers (update as published):
- Chapter 2 (germination/microbiome factors): *Scientific Reports* (in prep)
- Chapter 3 (full microbiome analysis): *Environmental Science & Technology* (in prep)

---

## Contact

Montoya Geoenvironmental Lab  
Department of Civil, Construction and Environmental Engineering  
NC State University  
PI: Dr. Brina Montoya
