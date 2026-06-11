# =============================================================================
# 00_config.R
# MICP Microbiome Analysis Pipeline — User Configuration
# =============================================================================
# PURPOSE:
#   This is the ONLY file you should need to edit between experiments.
#   All file paths, treatment labels, species codes, analysis parameters,
#   and output preferences are defined here and passed to all downstream
#   scripts automatically.
#
# USAGE:
#   Source this file at the top of any analysis script:
#     source("R/00_config.R")
#   or run the full pipeline via:
#     source("R/run_pipeline.R")
#
# NOTES FOR FUTURE USERS:
#   - Treatment naming must be consistent across OTU table, taxonomy table,
#     and metadata. If T-BAS or sequencing output changes column names,
#     update the column name mappings below (section: DATA COLUMN NAMES).
#   - The pipeline was developed for Nanopore full-length 16S data processed
#     through NanoClust/NanoScript → T-BAS (LIFE1 reference tree). If
#     switching to short-read (Illumina) data, rarefaction depth logic and
#     singleton filtering assumptions in notes may need adjustment.
#   - All outputs (figures, tables) are written to outputs/ automatically.
#     Set OVERWRITE_OUTPUTS = FALSE to protect existing results.
# =============================================================================


# -----------------------------------------------------------------------------
# EXPERIMENT METADATA
# -----------------------------------------------------------------------------
EXPERIMENT_ID   <- "ERDC_Task1"        # Used in output filenames & plot titles
EXPERIMENT_DESC <- "MICP Germination Phase — Task 1 Microbiome"
SEQUENCING_PLATFORM <- "Nanopore MinION (full-length 16S)"


# -----------------------------------------------------------------------------
# INPUT FILE PATHS
# -----------------------------------------------------------------------------
# Accepted formats: .csv for OTU/taxonomy/metadata; .nwk/.tre for tree.
# Use forward slashes or double backslashes on Windows.
# Relative paths are resolved from the repo root.

PATHS <- list(
  otu_table   = "data/raw/phyloseq_otu_table.csv",
  taxonomy    = "data/raw/phyloseq_taxonomy_table.csv",
  metadata    = "data/raw/ERDC_Metadata_corrected.csv",
  tree        = "data/raw/annotated_tree.nwk"
)

# Set to NULL to skip tree-based analyses (UniFrac, Faith's PD):
#   PATHS$tree <- NULL


# -----------------------------------------------------------------------------
# DATA COLUMN NAMES
# -----------------------------------------------------------------------------
# Map your actual CSV column names to the internal names used by the pipeline.
# Change the RIGHT side of each pair if your column headers differ.

COL_NAMES <- list(
  sample_id      = "Sample_Code",     # Row identifier in metadata & OTU table
  treatment      = "Treatment",        # Treatment group column in metadata
  species        = "Species",          # Plant species column in metadata
  otu_feature_id = "FeatureID",        # Row identifier in taxonomy table
  tax_domain     = "DOMAIN"            # Renamed internally to "Kingdom"
)


# -----------------------------------------------------------------------------
# TREATMENT LABELS
# -----------------------------------------------------------------------------
# Define the mapping from raw labels in your metadata CSV to the display
# labels used in all plots and statistical output.
# Order determines factor level order (left to right in plots).

TREATMENT_MAP <- list(
  raw_levels     = c("Control", "MICP 1", "MICP 2", "MICP 3"),
  display_labels = c("Control", "U250",   "U100+P", "U50+P")
)

# Color palette for treatment groups (must match length of display_labels)
# Default uses RColorBrewer "Set1"; override here with hex codes if needed:
#   TREATMENT_COLORS <- c("#000000", "#E41A1C", "#377EB8", "#4DAF4A")
TREATMENT_COLORS <- NULL   # NULL = use scale_color_brewer(palette = "Set1")


# -----------------------------------------------------------------------------
# SAMPLES TO EXCLUDE
# -----------------------------------------------------------------------------
# List sample IDs (using Sample_Code values) to drop before analysis.
# Common reasons: failed library prep, outlier contamination, incomplete data.
# Set to character(0) to include all samples.

EXCLUDE_SAMPLES <- c(
  "FEAR_C_5", "FEAR_M1_6", "FEAR_M2_7", "FEAR_M3_8"
  # Add future exclusions here, e.g.:
  # "TRAE_M2_11"
)


# -----------------------------------------------------------------------------
# SAMPLE ID CORRECTIONS
# -----------------------------------------------------------------------------
# Named vector of regex find → replacement pairs applied to Sample_Code
# before any analysis. Used to fix upstream typos from T-BAS or sequencing.
# Set to NULL to skip.

SAMPLE_ID_FIXES <- c(
  "^SEIT_" = "SIET_"   # T-BAS output uses SIET; metadata uses SEIT
  # Add future corrections here
)


# -----------------------------------------------------------------------------
# PLANT SPECIES — MICP PRE-ADAPTATION CLASSIFICATIONS
# -----------------------------------------------------------------------------
# Groups based on documented stress tolerance relevant to MICP field conditions:
#   - "Dual tolerant":      drought + alkalinity
#   - "Drought tolerant":   drought only
#   - "Alkalinity tolerant": alkalinity/high pH only
#   - "Neither":            no documented pre-adaptation to MICP stressors
#   - "No plants":          unplanted control pots
#
# Update this list when adding new species to the experiment.
# Any species code not listed here will receive NA for Adaptation_Group.

ADAPTATION_GROUPS <- c(
  DICL = "Dual tolerant", CYDA = "Dual tolerant", ERCU = "Dual tolerant",
  FEAR = "Dual tolerant", AGST = "Dual tolerant", COVA = "Dual tolerant",
  LECA = "Dual tolerant",
  PANO = "Drought tolerant", SECE = "Drought tolerant",
  SCSC = "Alkalinity tolerant", FERU = "Alkalinity tolerant",
  FEBR = "Alkalinity tolerant", POPR = "Alkalinity tolerant",
  CAVU = "Alkalinity tolerant", CHFA = "Alkalinity tolerant",
  SEIT = "Neither", CHLA = "Neither", ELVI = "Neither", TRAE = "Neither",
  NOPLANTS = "No plants"
)

# Shapes for adaptation groups in ordination plots
ADAPTATION_SHAPES <- c(
  "Dual tolerant"       = 16,   # filled circle
  "Drought tolerant"    = 17,   # filled triangle up
  "Alkalinity tolerant" = 15,   # filled square
  "Neither"             = 4,    # cross / X
  "No plants"           = 8     # asterisk
)


# -----------------------------------------------------------------------------
# ALPHA DIVERSITY SETTINGS
# -----------------------------------------------------------------------------
ALPHA_METRICS <- c("Observed", "Shannon", "Simpson")
# Faith's PD is always computed when a tree is available (PATHS$tree != NULL)

ALPHA_POSTHOC_METHOD <- "bonferroni"   # Correction for Dunn pairwise tests


# -----------------------------------------------------------------------------
# BETA DIVERSITY & RAREFACTION SETTINGS
# -----------------------------------------------------------------------------
RAREFY_SEED   <- 42       # Random seed for reproducibility
RAREFY_DEPTH  <- NULL     # NULL = auto (minimum sample depth); or set integer

BETA_DISTANCES    <- c("bray", "unifrac", "wunifrac")
PERMANOVA_PERMS   <- 999
PERMANOVA_FORMULA <- "~ Treatment + Species"   # RHS only; lhs = distance matrix


# -----------------------------------------------------------------------------
# RELATIVE ABUNDANCE BAR PLOT SETTINGS
# -----------------------------------------------------------------------------
TOP_N_GENERA <- 12    # Number of most abundant genera to display

# High-contrast palette (up to 20 colors). Adjust length to match TOP_N_GENERA.
# Colors are assigned in descending abundance order.
GENUS_PALETTE <- c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#EA5E5F", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FFD400", "#B15928", "#FFFAC8", "#800000", "#AAFFC3",
  "#808000", "#FFD8B1", "#000075", "#A9A9A9", "#000000"
)


# -----------------------------------------------------------------------------
# N-CYCLING GUILD DEFINITIONS
# -----------------------------------------------------------------------------
# Genera are assigned to guilds based on literature-documented N-transformation
# roles. Multi-guild membership is permitted (e.g., Paenibacillus contributes
# to Ureolytic, Denitrifier, AND Diazotroph guilds independently).
# Guild-level relative abundances are therefore NOT mutually exclusive.
#
# To add a new guild, append a new named entry to this list.
# To update membership, add/remove genus names within each vector.

N_GUILDS <- list(
  Ureolytic = c(
    "Sporosarcina", "Bacillus", "Paenibacillus", "Lysinibacillus",
    "Neobacillus", "Priestia", "Gottfriedia", "Lederbergia",
    "Ectobacillus", "Paenisporosarcina"
  ),
  AOB = c(
    "Nitrosomonas", "Nitrosospira", "Nitrosovibrio"
  ),
  NOB = c(
    "Nitrospira", "Nitrobacter", "Nitrotoga"
  ),
  Denitrifier = c(
    "Alcaligenes", "Castellaniella", "Delftia", "Hyphomicrobium",
    "Simplicispira", "Bacillus", "Paenibacillus", "Lysinibacillus",
    "Aminobacter", "Bosea", "Brevundimonas", "Caballeronia",
    "Enterobacter", "Massilia", "Mycoplana", "Pantoea",
    "Rhodopseudomonas", "Shinella", "Thauera", "Dechloromonas",
    "Aromatoleum", "Sulfurimonas", "Paracoccus"
  ),
  Diazotroph = c(
    "Mesorhizobium", "Noviherbaspirillum", "Paenibacillus",
    "Clostridium", "Microchaete", "Aminobacter", "Caballeronia",
    "Devosia", "Enterobacter", "Pantoea", "Rhodopseudomonas",
    "Rhizobium", "Bradyrhizobium", "Azospirillum",
    "Herbaspirillum", "Azoarcus", "Frankia"
  )
)

GUILD_POSTHOC_METHOD <- "BH"   # Benjamini-Hochberg for guild Dunn tests


# -----------------------------------------------------------------------------
# OUTPUT SETTINGS
# -----------------------------------------------------------------------------
OUT_DIR_FIGS   <- "outputs/figures"
OUT_DIR_TABLES <- "outputs/tables"

# Figure export formats. Options: "png", "tiff", "pdf", "svg"
# TIFF uses LZW compression (journal submission standard)
FIGURE_FORMATS <- c("tiff", "png")
FIGURE_DPI     <- 300

# Set FALSE to skip overwriting files that already exist in outputs/
OVERWRITE_OUTPUTS <- TRUE


# -----------------------------------------------------------------------------
# END CONFIG
# -----------------------------------------------------------------------------
cat(sprintf(
  "[config] Loaded: %s | Experiment: %s\n",
  format(Sys.time(), "%Y-%m-%d %H:%M"),
  EXPERIMENT_ID
))
