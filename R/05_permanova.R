# =============================================================================
# 05_permanova.R
# MICP Microbiome Analysis Pipeline — PERMANOVA & Betadisper
# =============================================================================
# PURPOSE:
#   Tests whether community composition differs significantly across Treatment
#   and Species using PERMANOVA (adonis2, by = "terms"). Also assesses whether
#   within-group dispersion (not just centroid location) differs using
#   betadisper/PERMDISP, which is an important validity check for PERMANOVA
#   interpretation.
#
# INPUT:  data/processed/ps_rarefied.rds
# OUTPUT: outputs/tables/permanova_results.csv
#         outputs/tables/pairwise_permanova_*.csv
#         outputs/figures/betadisper_boxplot.{tiff,png}
#
# DEPENDENCIES: tidyverse, phyloseq, vegan
#
# NOTES:
#   - adonis2(by = "terms") is used throughout. This tests predictors
#     sequentially (Treatment first, then Species), so results are
#     order-dependent. The formula is set in PERMANOVA_FORMULA (00_config.R).
#   - PERMDISP is a validity check for PERMANOVA, not a standalone result.
#     A significant PERMDISP indicates unequal dispersions, which can
#     inflate PERMANOVA significance — report accordingly.
#   - Pairwise PERMANOVA uses Bonferroni correction. Results with
#     p_bonferroni < 0.05 should be interpreted with the PERMDISP result
#     in mind.
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


# -----------------------------------------------------------------------------
# 1. LOAD & PREPARE METADATA
# -----------------------------------------------------------------------------
cat("[05] Loading rarefied phyloseq object...\n")
ps_rare <- readRDS("data/processed/ps_rarefied.rds")

has_tree <- !is.null(phy_tree(ps_rare, errorIfNULL = FALSE))

samp_df <- data.frame(sample_data(ps_rare))
samp_df$Treatment <- factor(samp_df$Treatment,
                             levels = TREATMENT_MAP$display_labels)
samp_df$Species   <- factor(samp_df$Species)

cat("Treatment levels:", paste(levels(samp_df$Treatment), collapse = ", "), "\n")
cat("Species levels:",   paste(levels(samp_df$Species),   collapse = ", "), "\n")
cat("N samples:", nrow(samp_df), "\n\n")


# -----------------------------------------------------------------------------
# 2. COMPUTE DISTANCE MATRICES
# -----------------------------------------------------------------------------
cat("[05] Computing distance matrices...\n")

dist_list <- list(
  bray = phyloseq::distance(ps_rare, method = "bray")
)

if (has_tree) {
  dist_list$unifrac  <- suppressWarnings(
    phyloseq::distance(ps_rare, method = "unifrac")
  )
  dist_list$wunifrac <- suppressWarnings(
    phyloseq::distance(ps_rare, method = "wunifrac")
  )
} else {
  cat("[05] No tree — UniFrac distances skipped.\n")
}

# Only run distances that were configured
dist_list <- dist_list[intersect(names(dist_list), BETA_DISTANCES)]


# -----------------------------------------------------------------------------
# 3. PERMANOVA (adonis2)
# -----------------------------------------------------------------------------
# Formula RHS comes from config; LHS is the distance matrix
formula_rhs <- PERMANOVA_FORMULA   # e.g. "~ Treatment + Species"
permanova_results <- list()

cat(sprintf("[05] Running PERMANOVA (formula: dist %s, perms = %d)\n",
            formula_rhs, PERMANOVA_PERMS))

for (dist_name in names(dist_list)) {
  cat(sprintf("\n=== PERMANOVA (%s) ===\n", dist_name))
  formula_full <- as.formula(paste("dist_list[[dist_name]]", formula_rhs))
  res <- adonis2(formula_full, data = samp_df,
                 permutations = PERMANOVA_PERMS, by = "terms")
  print(res)
  permanova_results[[dist_name]] <- res
}

# Export as flat CSV
perm_export <- map_dfr(names(permanova_results), function(d) {
  r <- permanova_results[[d]]
  tibble(
    distance  = d,
    term      = rownames(r),
    df        = r$Df,
    SumOfSqs  = round(r$SumOfSqs, 4),
    R2        = round(r$R2, 4),
    F         = round(r$F, 4),
    p         = r$`Pr(>F)`
  )
})
write.csv(perm_export,
          file.path(OUT_DIR_TABLES, "permanova_results.csv"),
          row.names = FALSE)
cat(sprintf("\n[05] PERMANOVA results saved to %s/permanova_results.csv\n",
            OUT_DIR_TABLES))


# -----------------------------------------------------------------------------
# 4. PAIRWISE PERMANOVA (Treatment contrasts with Bonferroni correction)
# -----------------------------------------------------------------------------
pairwise_permanova <- function(dist_matrix, grouping,
                               permutations = PERMANOVA_PERMS) {
  groups  <- levels(factor(grouping))
  pairs   <- combn(groups, 2, simplify = FALSE)

  map_dfr(pairs, function(pair) {
    keep     <- grouping %in% pair
    sub_grp  <- grouping[keep]
    sub_dist <- as.dist(as.matrix(dist_matrix)[keep, keep])
    res      <- adonis2(sub_dist ~ sub_grp, permutations = permutations,
                        by = "terms")
    tibble(
      Group1 = pair[1], Group2 = pair[2],
      F      = round(res$F[1], 4),
      R2     = round(res$R2[1], 4),
      p      = res$`Pr(>F)`[1]
    )
  }) %>%
    mutate(
      p_bonferroni = p.adjust(p, method = "bonferroni"),
      sig = case_when(
        p_bonferroni < 0.001 ~ "***",
        p_bonferroni < 0.01  ~ "**",
        p_bonferroni < 0.05  ~ "*",
        p_bonferroni < 0.1   ~ ".",
        TRUE                 ~ "ns"
      )
    )
}

for (dist_name in names(dist_list)) {
  cat(sprintf("\n=== Pairwise PERMANOVA (%s) ===\n", dist_name))
  pw <- pairwise_permanova(dist_list[[dist_name]], samp_df$Treatment)
  print(pw)
  write.csv(pw,
            file.path(OUT_DIR_TABLES,
                      sprintf("pairwise_permanova_%s.csv", dist_name)),
            row.names = FALSE)
}


# -----------------------------------------------------------------------------
# 5. BETADISPER (PERMDISP)
# -----------------------------------------------------------------------------
# NOTE: PERMDISP is reported as a validity check for PERMANOVA, not as a
# primary result. A significant result means dispersions are unequal, which
# should be acknowledged when interpreting PERMANOVA significance.

cat("\n=== Betadisper: PERMDISP (Bray-Curtis, Treatment) ===\n")
bd_treat <- betadisper(dist_list$bray, samp_df$Treatment)
print(permutest(bd_treat))

cat("\n=== Betadisper: PERMDISP (Bray-Curtis, Species) ===\n")
bd_spec <- betadisper(dist_list$bray, samp_df$Species)
print(permutest(bd_spec))

# Boxplot of distances-to-centroid by Treatment
bd_df <- data.frame(
  Distance  = bd_treat$distances,
  Treatment = samp_df$Treatment[match(names(bd_treat$distances),
                                      rownames(samp_df))]
)

p_bd <- ggplot(bd_df, aes(x = Treatment, y = Distance, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Distance to Centroid by Treatment (Bray-Curtis)",
    y     = "Distance to centroid",
    x     = "Treatment",
    caption = "Betadisper/PERMDISP — validity check for PERMANOVA"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(p_bd)
save_figure(p_bd, "betadisper_boxplot", width = 7, height = 5)

cat("[05] PERMANOVA analysis complete.\n")
