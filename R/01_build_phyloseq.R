# =============================================================================
# 01_build_phyloseq.R
# MICP Microbiome Analysis Pipeline — Data Import & phyloseq Construction
# =============================================================================
# PURPOSE:
#   Loads raw OTU table, taxonomy table, sample metadata, and phylogenetic
#   tree; applies sample exclusions and ID corrections defined in config;
#   builds and validates the phyloseq object used by all downstream scripts.
#
# INPUT:  Paths defined in PATHS (00_config.R)
# OUTPUT: ps.rds saved to data/processed/
#
# DEPENDENCIES: tidyverse, phyloseq, ape
# =============================================================================

source("R/00_config.R")

library(tidyverse)
library(phyloseq)


# -----------------------------------------------------------------------------
# HELPER: safe file reader with informative errors
# -----------------------------------------------------------------------------
read_input <- function(path, type = "csv") {
  if (is.null(path)) return(NULL)
  if (!file.exists(path)) {
    stop(sprintf(
      "[01_build_phyloseq] File not found: %s\n  Check PATHS in 00_config.R",
      path
    ))
  }
  if (type == "csv") read_csv(path, show_col_types = FALSE)
  else if (type == "tree") ape::read.tree(path)
}


# -----------------------------------------------------------------------------
# 1. LOAD RAW DATA
# -----------------------------------------------------------------------------
cat("[01] Loading input files...\n")

otu_raw  <- read_input(PATHS$otu_table)
tax_raw  <- read_input(PATHS$taxonomy)
meta_raw <- read_input(PATHS$metadata)
tree     <- read_input(PATHS$tree, type = "tree")

if (is.null(tree)) {
  cat("[01] No tree file provided — UniFrac and Faith's PD will be skipped.\n")
}


# -----------------------------------------------------------------------------
# 2. OTU TABLE
# -----------------------------------------------------------------------------
# Expected format from T-BAS: samples as rows, OTUs as columns, first column
# is the sample identifier (COL_NAMES$sample_id).
# We transpose to taxa-as-rows for phyloseq convention.

otu_mat <- otu_raw %>%
  column_to_rownames(COL_NAMES$sample_id) %>%
  as.matrix() %>%
  t()   # phyloseq: rows = OTUs, columns = samples

cat(sprintf("[01] OTU table: %d OTUs × %d samples\n",
            nrow(otu_mat), ncol(otu_mat)))


# -----------------------------------------------------------------------------
# 3. TAXONOMY TABLE
# -----------------------------------------------------------------------------
# Expected format: FeatureID column + one column per taxonomic rank.
# The DOMAIN column is renamed to Kingdom for phyloseq compatibility.

tax_mat <- tax_raw %>%
  column_to_rownames(COL_NAMES$otu_feature_id) %>%
  as.matrix()

# Rename DOMAIN → Kingdom if present
if (COL_NAMES$tax_domain %in% colnames(tax_mat)) {
  colnames(tax_mat)[colnames(tax_mat) == COL_NAMES$tax_domain] <- "Kingdom"
}

cat(sprintf("[01] Taxonomy table: %d features, ranks: %s\n",
            nrow(tax_mat), paste(colnames(tax_mat), collapse = ", ")))


# -----------------------------------------------------------------------------
# 4. SAMPLE METADATA
# -----------------------------------------------------------------------------
# Apply sample exclusions and ID corrections before building phyloseq.

meta_clean <- meta_raw

# Apply sample ID corrections (regex find → replace on sample ID column)
if (!is.null(SAMPLE_ID_FIXES)) {
  for (pattern in names(SAMPLE_ID_FIXES)) {
    meta_clean[[COL_NAMES$sample_id]] <- gsub(
      pattern,
      SAMPLE_ID_FIXES[[pattern]],
      meta_clean[[COL_NAMES$sample_id]]
    )
  }
  cat(sprintf("[01] Applied %d sample ID correction(s).\n",
              length(SAMPLE_ID_FIXES)))
}

# Drop excluded samples
n_before <- nrow(meta_clean)
meta_clean <- meta_clean %>%
  filter(!.data[[COL_NAMES$sample_id]] %in% EXCLUDE_SAMPLES)
n_dropped <- n_before - nrow(meta_clean)
if (n_dropped > 0) {
  cat(sprintf("[01] Excluded %d sample(s): %s\n",
              n_dropped, paste(EXCLUDE_SAMPLES, collapse = ", ")))
}

# Remap treatment labels (raw → display)
meta_clean[[COL_NAMES$treatment]] <- factor(
  meta_clean[[COL_NAMES$treatment]],
  levels = TREATMENT_MAP$raw_levels,
  labels = TREATMENT_MAP$display_labels
)

meta_clean <- meta_clean %>%
  column_to_rownames(COL_NAMES$sample_id)

cat(sprintf("[01] Metadata: %d samples after filtering.\n", nrow(meta_clean)))


# -----------------------------------------------------------------------------
# 5. ALIGN SAMPLES ACROSS ALL COMPONENTS
# -----------------------------------------------------------------------------
# Ensure sample IDs match between OTU table and metadata.
# Mismatches here are the most common source of phyloseq build failures.

shared_samples <- intersect(colnames(otu_mat), rownames(meta_clean))
only_otu  <- setdiff(colnames(otu_mat), rownames(meta_clean))
only_meta <- setdiff(rownames(meta_clean), colnames(otu_mat))

if (length(only_otu) > 0) {
  warning(sprintf(
    "[01] %d sample(s) in OTU table but NOT in metadata — will be dropped:\n  %s",
    length(only_otu), paste(only_otu, collapse = ", ")
  ))
}
if (length(only_meta) > 0) {
  warning(sprintf(
    "[01] %d sample(s) in metadata but NOT in OTU table — will be dropped:\n  %s",
    length(only_meta), paste(only_meta, collapse = ", ")
  ))
}

otu_mat    <- otu_mat[, shared_samples, drop = FALSE]
meta_clean <- meta_clean[shared_samples, , drop = FALSE]

cat(sprintf("[01] Samples aligned: %d shared between OTU table and metadata.\n",
            length(shared_samples)))


# -----------------------------------------------------------------------------
# 6. BUILD PHYLOSEQ OBJECT
# -----------------------------------------------------------------------------
OTU  <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX  <- tax_table(tax_mat)
SAMP <- sample_data(meta_clean)

if (!is.null(tree)) {
  ps <- phyloseq(OTU, TAX, SAMP, tree)
} else {
  ps <- phyloseq(OTU, TAX, SAMP)
}

# Remove OTUs with zero total counts (can arise after sample exclusions)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Add adaptation group classifications
sample_data(ps)$Adaptation_Group <- ADAPTATION_GROUPS[
  as.character(sample_data(ps)[[COL_NAMES$species]])
]


# -----------------------------------------------------------------------------
# 7. VALIDATION REPORT
# -----------------------------------------------------------------------------
cat("\n========== phyloseq Object Summary ==========\n")
cat(sprintf("  Samples : %d\n", nsamples(ps)))
cat(sprintf("  OTUs    : %d (after pruning zero-count taxa)\n", ntaxa(ps)))
cat(sprintf("  Tree    : %s\n", ifelse(!is.null(tree), "present", "absent")))
cat(sprintf("  Ranks   : %s\n\n", paste(rank_names(ps), collapse = ", ")))

cat("Samples per treatment:\n")
print(table(sample_data(ps)[[COL_NAMES$treatment]]))

cat("\nSamples per species:\n")
print(table(sample_data(ps)[[COL_NAMES$species]]))

cat("\nSample depth range:\n")
depths <- sample_sums(ps)
cat(sprintf("  Min: %d | Median: %.0f | Max: %d\n\n",
            min(depths), median(depths), max(depths)))

# Warn if any treatment level is missing (factor levels with n=0 cause issues)
missing_treatments <- names(which(table(sample_data(ps)[[COL_NAMES$treatment]]) == 0))
if (length(missing_treatments) > 0) {
  warning(sprintf(
    "[01] Treatment(s) with NO samples after filtering: %s\n  Check EXCLUDE_SAMPLES and TREATMENT_MAP in 00_config.R",
    paste(missing_treatments, collapse = ", ")
  ))
}


# -----------------------------------------------------------------------------
# 8. SAVE PROCESSED OBJECT
# -----------------------------------------------------------------------------
out_path <- "data/processed/ps.rds"
saveRDS(ps, out_path)
cat(sprintf("[01] phyloseq object saved to: %s\n", out_path))
