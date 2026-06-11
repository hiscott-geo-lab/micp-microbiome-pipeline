# =============================================================================
# 02_relative_abundance.R
# MICP Microbiome Analysis Pipeline — Stacked Bar Plots (Relative Abundance)
# =============================================================================
# PURPOSE:
#   Generates genus-level relative abundance bar plots:
#     A) Averaged by Treatment (one bar per treatment group)
#     B) Faceted by Species × Treatment (one panel per plant species)
#
# INPUT:  data/processed/ps.rds
# OUTPUT: outputs/figures/bar_treatment.{tiff,png}
#         outputs/figures/bar_by_species.{tiff,png}
#
# DEPENDENCIES: tidyverse, phyloseq, scales
# =============================================================================

source("R/00_config.R")

library(tidyverse)
library(phyloseq)
library(scales)


# -----------------------------------------------------------------------------
# HELPER: save figure in all configured formats
# -----------------------------------------------------------------------------
save_figure <- function(plot, filename_stem, width, height) {
  for (fmt in FIGURE_FORMATS) {
    out_path <- file.path(OUT_DIR_FIGS, paste0(filename_stem, ".", fmt))
    if (!OVERWRITE_OUTPUTS && file.exists(out_path)) {
      cat(sprintf("[save] Skipping (exists): %s\n", out_path))
      next
    }
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
# 1. LOAD & PREPARE
# -----------------------------------------------------------------------------
cat("[02] Loading phyloseq object...\n")
ps <- readRDS("data/processed/ps.rds")

# Agglomerate to genus and convert to relative abundance
ps_genus <- tax_glom(ps, taxrank = "GENUS", NArm = FALSE)
ps_rel   <- transform_sample_counts(ps_genus, function(x) x / sum(x))

# Select top N genera by mean relative abundance across all samples
top_genera <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:TOP_N_GENERA]
ps_top     <- prune_taxa(top_genera, ps_rel)

# Melt to long format for ggplot
bar_df <- psmelt(ps_top) %>%
  mutate(
    GENUS     = replace_na(GENUS, "Unclassified"),
    Treatment = factor(Treatment,
                       levels = TREATMENT_MAP$display_labels)
  )

# Build genus order (descending total abundance) for consistent palette assignment
genus_order <- bar_df %>%
  group_by(GENUS) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(GENUS)

# Assign palette — pad with grey if more genera than palette colors
n_genera <- length(genus_order)
palette_use <- GENUS_PALETTE
if (n_genera > length(palette_use)) {
  palette_use <- c(palette_use, rep("#CCCCCC", n_genera - length(palette_use)))
  warning(sprintf(
    "[02] More genera (%d) than palette colors (%d). Extra genera shown in grey.",
    n_genera, length(GENUS_PALETTE)
  ))
}
names(palette_use) <- genus_order

bar_df <- bar_df %>%
  mutate(GENUS = factor(GENUS, levels = genus_order))

cat(sprintf("[02] Top %d genera selected. Most abundant: %s\n",
            TOP_N_GENERA, paste(genus_order[1:min(5, n_genera)], collapse = ", ")))


# -----------------------------------------------------------------------------
# 2. PLOT A — Averaged by Treatment
# -----------------------------------------------------------------------------
bar_treatment <- bar_df %>%
  group_by(Treatment, GENUS) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

p_bar_treatment <- ggplot(bar_treatment,
                          aes(x = Treatment, y = Abundance, fill = GENUS)) +
  geom_bar(stat = "identity", position = "fill",
           color = "white", linewidth = 0.3) +
  scale_fill_manual(values = palette_use) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = sprintf("Relative Abundance by Treatment — Top %d Genera", TOP_N_GENERA),
    x     = "Treatment",
    y     = "Mean Relative Abundance",
    fill  = "Genus"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(hjust = 1),
    legend.text  = element_text(face = "italic")
  )

print(p_bar_treatment)
save_figure(p_bar_treatment, "bar_treatment", width = 10, height = 6)


# -----------------------------------------------------------------------------
# 3. PLOT B — Faceted by Species
# -----------------------------------------------------------------------------
# ncol adjusts automatically: 5 columns for ≤20 species; fewer if needed
n_species <- length(unique(bar_df$Species))
facet_ncol <- min(5, n_species)

p_bar_species <- ggplot(bar_df,
                        aes(x = Treatment, y = Abundance, fill = GENUS)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = palette_use) +
  facet_wrap(~ Species, ncol = facet_ncol) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = sprintf("Relative Abundance by Species × Treatment — Top %d Genera",
                    TOP_N_GENERA),
    x     = "Treatment",
    y     = "Relative Abundance",
    fill  = "Genus"
  ) +
  theme_bw() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    strip.text      = element_text(size = 8),
    legend.position = "bottom",
    legend.text     = element_text(face = "italic")
  )

print(p_bar_species)
save_figure(p_bar_species, "bar_by_species", width = 22, height = 12)

cat("[02] Relative abundance plots complete.\n")
