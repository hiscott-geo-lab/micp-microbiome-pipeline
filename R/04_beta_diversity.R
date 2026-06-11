# =============================================================================
# 04_beta_diversity.R
# MICP Microbiome Analysis Pipeline — Beta Diversity Ordination
# =============================================================================
# PURPOSE:
#   Computes and visualizes beta diversity using:
#     - Bray-Curtis dissimilarity (PCoA and NMDS)
#     - Unweighted UniFrac (PCoA)  [requires tree]
#     - Weighted UniFrac (PCoA)    [requires tree]
#   Points are colored by Treatment and shaped by plant tolerance classification
#   (Adaptation_Group) to facilitate visualization of MICP pre-adaptation.
#
# INPUT:  data/processed/ps.rds
# OUTPUT: outputs/figures/beta_*.{tiff,png}
#         Rarefied phyloseq saved to data/processed/ps_rarefied.rds
#
# DEPENDENCIES: tidyverse, phyloseq, vegan
#
# NOTES ON RAREFACTION:
#   Rarefaction to the minimum sample depth is applied before ordination to
#   control for sequencing depth differences. The seed is set in 00_config.R
#   for reproducibility. If RAREFY_DEPTH is set to NULL in config, the
#   minimum observed depth is used. Setting an explicit depth is recommended
#   when one or more samples have unusually low read counts that would
#   otherwise force all samples to a very low rarefaction cutoff.
# =============================================================================

source("R/00_config.R")

library(tidyverse)
library(phyloseq)
library(vegan)


# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
save_figure <- function(plot, filename_stem, width, height) {
  for (fmt in FIGURE_FORMATS) {
    out_path <- file.path(OUT_DIR_FIGS, paste0(filename_stem, ".", fmt))
    if (!OVERWRITE_OUTPUTS && file.exists(out_path)) next
    if (fmt == "tiff") {
      ggsave(out_path, plot, width = width, height = height,
             dpi = FIGURE_DPI, compression = "lzw")
    } else {
      ggsave(out_path, plot, width = width, height = height, dpi = FIGURE_DPI)
    }
    cat(sprintf("[save] Written: %s\n", out_path))
  }
}

ordination_theme <- function() {
  list(
    scale_color_brewer(palette = "Set1"),
    scale_shape_manual(values = ADAPTATION_SHAPES, drop = FALSE),
    theme_bw(),
    labs(color = "Treatment", shape = "Plant Tolerance Traits")
  )
}


# -----------------------------------------------------------------------------
# 1. LOAD & RAREFY
# -----------------------------------------------------------------------------
cat("[04] Loading phyloseq object...\n")
ps <- readRDS("data/processed/ps.rds")

has_tree <- !is.null(phy_tree(ps, errorIfNULL = FALSE))

depths     <- sample_sums(ps)
rare_depth <- RAREFY_DEPTH %||% min(depths)
`%||%`     <- function(a, b) if (!is.null(a)) a else b
rare_depth <- if (!is.null(RAREFY_DEPTH)) RAREFY_DEPTH else min(depths)

cat(sprintf("[04] Rarefying to depth: %d (seed = %d)\n", rare_depth, RAREFY_SEED))
cat(sprintf("[04] Samples below rarefy depth (will be dropped): %d\n",
            sum(depths < rare_depth)))

ps_rare <- rarefy_even_depth(ps,
                             sample.size = rare_depth,
                             rngseed     = RAREFY_SEED,
                             replace     = FALSE,
                             verbose     = FALSE)

# Re-attach adaptation groups (rarefy drops sample_data slots sometimes)
sample_data(ps_rare)$Adaptation_Group <- ADAPTATION_GROUPS[
  as.character(sample_data(ps_rare)[[COL_NAMES$species]])
]

saveRDS(ps_rare, "data/processed/ps_rarefied.rds")
cat(sprintf("[04] Rarefied ps saved to data/processed/ps_rarefied.rds (%d samples)\n",
            nsamples(ps_rare)))


# -----------------------------------------------------------------------------
# 2. BRAY-CURTIS PCoA
# -----------------------------------------------------------------------------
cat("[04] Computing Bray-Curtis PCoA...\n")
ord_bc_pcoa <- ordinate(ps_rare, method = "PCoA", distance = "bray")

p_bc_pcoa <- plot_ordination(ps_rare, ord_bc_pcoa, color = "Treatment") +
  geom_point(aes(shape = Adaptation_Group), size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Treatment), linetype = 2) +
  ordination_theme() +
  labs(title = "Bray-Curtis PCoA")

print(p_bc_pcoa)
save_figure(p_bc_pcoa, "beta_braycurtis_pcoa", width = 10, height = 7)


# -----------------------------------------------------------------------------
# 3. BRAY-CURTIS NMDS
# -----------------------------------------------------------------------------
cat("[04] Computing Bray-Curtis NMDS (trymax = 100)...\n")
ord_bc_nmds <- ordinate(ps_rare, method = "NMDS", distance = "bray",
                        trymax = 100)

p_bc_nmds <- plot_ordination(ps_rare, ord_bc_nmds, color = "Treatment") +
  geom_point(aes(shape = Adaptation_Group), size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Treatment), linetype = 2) +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("Stress = %.3f", ord_bc_nmds$stress),
           hjust = 1.1, vjust = -0.5, size = 3.5) +
  ordination_theme() +
  labs(title = "Bray-Curtis NMDS")

print(p_bc_nmds)
save_figure(p_bc_nmds, "beta_braycurtis_nmds", width = 10, height = 7)

# NMDS stress diagnostic note
if (ord_bc_nmds$stress > 0.2) {
  warning(sprintf(
    "[04] NMDS stress = %.3f (> 0.2). Ordination may not reliably represent distances. Consider increasing trymax or inspecting outlier samples.",
    ord_bc_nmds$stress
  ))
}


# -----------------------------------------------------------------------------
# 4. UniFrac PCoA (requires tree)
# -----------------------------------------------------------------------------
if (!has_tree) {
  cat("[04] Skipping UniFrac ordinations (no tree).\n")
} else {
  # --- Unweighted UniFrac ---
  cat("[04] Computing Unweighted UniFrac PCoA...\n")
  ord_uu <- suppressWarnings(
    ordinate(ps_rare, method = "PCoA", distance = "unifrac")
  )

  p_uu <- plot_ordination(ps_rare, ord_uu, color = "Treatment") +
    geom_point(aes(shape = Adaptation_Group), size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = Treatment), linetype = 2) +
    ordination_theme() +
    labs(title = "Unweighted UniFrac PCoA")

  print(p_uu)
  save_figure(p_uu, "beta_unifrac_unweighted_pcoa", width = 10, height = 7)

  # --- Weighted UniFrac ---
  cat("[04] Computing Weighted UniFrac PCoA...\n")
  ord_wu <- suppressWarnings(
    ordinate(ps_rare, method = "PCoA", distance = "wunifrac")
  )

  p_wu <- plot_ordination(ps_rare, ord_wu, color = "Treatment") +
    geom_point(aes(shape = Adaptation_Group), size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = Treatment), linetype = 2) +
    ordination_theme() +
    labs(title = "Weighted UniFrac PCoA")

  print(p_wu)
  save_figure(p_wu, "beta_unifrac_weighted_pcoa", width = 10, height = 7)
}

cat("[04] Beta diversity ordination complete.\n")
