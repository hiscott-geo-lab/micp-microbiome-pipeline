# =============================================================================
# 03_alpha_diversity.R
# MICP Microbiome Analysis Pipeline — Alpha Diversity
# =============================================================================
# PURPOSE:
#   Calculates and visualizes alpha diversity for Shannon, Simpson, Observed
#   OTUs, and Faith's Phylogenetic Diversity (when a tree is available).
#   Runs Kruskal-Wallis omnibus tests and Dunn post-hoc pairwise comparisons.
#
# INPUT:  data/processed/ps.rds
# OUTPUT: outputs/figures/alpha_diversity.{tiff,png}
#         outputs/figures/alpha_faiths_pd.{tiff,png}
#         outputs/figures/alpha_faiths_pd_by_species.{tiff,png}
#         outputs/tables/alpha_diversity_stats.csv
#         outputs/tables/alpha_faiths_pd_stats.csv
#
# DEPENDENCIES: tidyverse, phyloseq, picante, dunn.test
#
# NOTES:
#   Singletons are removed during the T-BAS pipeline (expected for Nanopore
#   consensus OTUs). Shannon and Simpson are robust to this; Observed OTU
#   counts will be conservative and should be interpreted accordingly.
#   Faith's PD requires a phylogenetic tree in the ps object.
# =============================================================================

source("R/00_config.R")

library(tidyverse)
library(phyloseq)


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

alpha_boxplot <- function(df, y_var, fill_var = "Treatment",
                          facet_var = NULL, ncol_facet = 5,
                          title = NULL, y_lab = NULL) {
  p <- ggplot(df, aes(x = .data[[fill_var]], y = .data[[y_var]],
                      fill = .data[[fill_var]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = title, x = fill_var, y = y_lab %||% y_var) +
    theme_bw() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(reformulate(facet_var), ncol = ncol_facet,
                        scales = "free_y")
  }
  p
}

# NULL coalescing (base R doesn't have %||%)
`%||%` <- function(a, b) if (!is.null(a)) a else b


# -----------------------------------------------------------------------------
# 1. LOAD
# -----------------------------------------------------------------------------
cat("[03] Loading phyloseq object...\n")
ps <- readRDS("data/processed/ps.rds")

has_tree <- !is.null(phy_tree(ps, errorIfNULL = FALSE))
cat(sprintf("[03] Phylogenetic tree: %s\n",
            ifelse(has_tree, "present — Faith's PD will be computed",
                   "absent  — Faith's PD will be skipped")))


# -----------------------------------------------------------------------------
# 2. SHANNON / SIMPSON / OBSERVED
# -----------------------------------------------------------------------------
cat("[03] Computing standard alpha diversity metrics...\n")

alpha_df <- suppressWarnings(
  estimate_richness(ps, measures = ALPHA_METRICS)
) %>%
  rownames_to_column("Sample_Code") %>%
  left_join(
    data.frame(sample_data(ps)) %>% rownames_to_column("Sample_Code"),
    by = "Sample_Code"
  )

# Kruskal-Wallis omnibus test
cat("\n=== Alpha Diversity: Kruskal-Wallis Tests ===\n")
kw_results <- map_dfr(ALPHA_METRICS, function(metric) {
  kw <- kruskal.test(reformulate(COL_NAMES$treatment, response = metric),
                     data = alpha_df)
  tibble(Metric = metric, chi_sq = round(kw$statistic, 3),
         df = kw$parameter, p = round(kw$p.value, 4))
})
print(kw_results)

# Dunn post-hoc pairwise comparisons
library(dunn.test)
cat(sprintf("\n=== Pairwise Dunn Tests (%s correction) ===\n",
            ALPHA_POSTHOC_METHOD))
dunn_results <- list()
for (metric in ALPHA_METRICS) {
  cat(sprintf("\n--- %s ---\n", metric))
  dunn_results[[metric]] <- dunn.test(
    alpha_df[[metric]],
    g      = alpha_df[[COL_NAMES$treatment]],
    method = ALPHA_POSTHOC_METHOD,
    kw     = FALSE,
    label  = TRUE,
    wrap   = TRUE
  )
}

# Save stats table
alpha_stats_out <- file.path(OUT_DIR_TABLES, "alpha_diversity_stats.csv")
write.csv(kw_results, alpha_stats_out, row.names = FALSE)
cat(sprintf("\n[03] KW results saved to: %s\n", alpha_stats_out))

# Plot: all standard metrics faceted
alpha_long <- alpha_df %>%
  pivot_longer(cols = all_of(ALPHA_METRICS),
               names_to = "Metric", values_to = "Value")

p_alpha <- ggplot(alpha_long, aes(x = Treatment, y = Value, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Alpha Diversity by Treatment",
    x     = "Treatment",
    y     = "Diversity Value"
  ) +
  theme_bw() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p_alpha)
save_figure(p_alpha, "alpha_diversity", width = 10, height = 5)


# -----------------------------------------------------------------------------
# 3. FAITH'S PHYLOGENETIC DIVERSITY
# -----------------------------------------------------------------------------
if (!has_tree) {
  cat("[03] Skipping Faith's PD (no tree).\n")
} else {
  library(picante)
  cat("[03] Computing Faith's Phylogenetic Diversity...\n")

  otu_for_pd  <- t(otu_table(ps))           # picante: samples × OTUs
  tree_for_pd <- phy_tree(ps)
  tree_for_pd <- prune.sample(otu_for_pd, tree_for_pd)

  pd_result <- pd(otu_for_pd, tree_for_pd, include.root = TRUE)

  pd_df <- pd_result %>%
    rownames_to_column("Sample_Code") %>%
    left_join(
      data.frame(sample_data(ps)) %>% rownames_to_column("Sample_Code"),
      by = "Sample_Code"
    )

  # Summary
  cat("\n=== Faith's PD Summary by Treatment ===\n")
  pd_summary <- pd_df %>%
    group_by(Treatment) %>%
    summarise(n = n(), mean_PD = round(mean(PD), 3), sd_PD = round(sd(PD), 3),
              min_PD = round(min(PD), 3), max_PD = round(max(PD), 3),
              .groups = "drop")
  print(pd_summary)
  write.csv(pd_summary,
            file.path(OUT_DIR_TABLES, "alpha_faiths_pd_stats.csv"),
            row.names = FALSE)

  # KW test
  kw_pd <- kruskal.test(PD ~ Treatment, data = pd_df)
  cat(sprintf(
    "\nFaith's PD — Kruskal-Wallis: chi-sq = %.3f, df = %d, p = %.4f\n",
    kw_pd$statistic, kw_pd$parameter, kw_pd$p.value
  ))

  # Dunn post-hoc
  cat(sprintf("\n--- Faith's PD (Dunn, %s) ---\n", ALPHA_POSTHOC_METHOD))
  dunn.test(pd_df$PD, g = pd_df$Treatment,
            method = ALPHA_POSTHOC_METHOD, kw = FALSE, label = TRUE, wrap = TRUE)

  # Plot: by Treatment
  p_pd <- alpha_boxplot(pd_df, y_var = "PD",
                        title = "Faith's Phylogenetic Diversity by Treatment",
                        y_lab = "Faith's PD")
  print(p_pd)
  save_figure(p_pd, "alpha_faiths_pd", width = 7, height = 5)

  # Plot: faceted by Species
  n_species  <- length(unique(pd_df$Species))
  facet_ncol <- min(5, n_species)

  p_pd_species <- alpha_boxplot(
    pd_df, y_var = "PD",
    facet_var  = "Species", ncol_facet = facet_ncol,
    title = "Faith's PD by Species × Treatment",
    y_lab = "Faith's PD"
  )
  print(p_pd_species)
  save_figure(p_pd_species, "alpha_faiths_pd_by_species", width = 22, height = 12)
}

cat("[03] Alpha diversity analysis complete.\n")
