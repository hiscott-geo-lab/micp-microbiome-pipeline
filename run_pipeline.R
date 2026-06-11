# =============================================================================
# run_pipeline.R
# MICP Microbiome Analysis Pipeline — Master Runner
# =============================================================================
# PURPOSE:
#   Runs the full pipeline end-to-end in the correct order.
#   This is the recommended entry point for new experiments.
#
# USAGE (from RStudio or R console):
#   setwd("/path/to/micp-microbiome-pipeline")
#   source("run_pipeline.R")
#
# TO RUN INDIVIDUAL STEPS ONLY:
#   source("R/00_config.R")      # Load config
#   source("R/01_build_phyloseq.R")  # Just the data import step
#
# PIPELINE ORDER:
#   00_config.R          — Configuration (sourced automatically by each script)
#   01_build_phyloseq.R  — Import data, build & save ps object
#   02_relative_abundance.R  — Stacked bar plots
#   03_alpha_diversity.R     — Alpha diversity metrics & stats
#   04_beta_diversity.R      — Ordination (PCoA, NMDS)
#   05_permanova.R           — PERMANOVA, pairwise PERMANOVA, PERMDISP
#   06_ncycling_guilds.R     — N-cycling functional guild analysis
#
# OUTPUTS:
#   data/processed/ps.rds           — Built phyloseq object
#   data/processed/ps_rarefied.rds  — Rarefied phyloseq object
#   outputs/figures/                — All plots (.tiff and .png)
#   outputs/tables/                 — All statistical result tables (.csv)
# =============================================================================

cat("=============================================================\n")
cat(" MICP Microbiome Pipeline\n")
cat(sprintf(" Started: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=============================================================\n\n")

# Ensure working directory is the repo root
if (!file.exists("R/00_config.R")) {
  stop("run_pipeline.R must be sourced from the repo root directory.\n",
       "  Use: setwd('/path/to/micp-microbiome-pipeline') then source('run_pipeline.R')")
}

# Create output directories if they don't exist
source("R/00_config.R")  # loads OUT_DIR_FIGS and OUT_DIR_TABLES
dir.create(OUT_DIR_FIGS,   showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_TABLES, showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# Run each step in order
steps <- list(
  "01 — Build phyloseq"         = "R/01_build_phyloseq.R",
  "02 — Relative abundance"     = "R/02_relative_abundance.R",
  "03 — Alpha diversity"        = "R/03_alpha_diversity.R",
  "04 — Beta diversity"         = "R/04_beta_diversity.R",
  "05 — PERMANOVA"              = "R/05_permanova.R",
  "06 — N-cycling guilds"       = "R/06_ncycling_guilds.R"
)

for (step_name in names(steps)) {
  cat(sprintf("\n----- %s -----\n", step_name))
  tryCatch(
    source(steps[[step_name]]),
    error = function(e) {
      cat(sprintf("\n[ERROR] Step failed: %s\n  %s\n", step_name, e$message))
      cat("  Pipeline halted. Fix the issue above, then re-run from this step.\n")
      stop(e)
    }
  )
}

cat("\n=============================================================\n")
cat(sprintf(" Pipeline complete: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf(" Figures : %s\n", OUT_DIR_FIGS))
cat(sprintf(" Tables  : %s\n", OUT_DIR_TABLES))
cat("=============================================================\n")
