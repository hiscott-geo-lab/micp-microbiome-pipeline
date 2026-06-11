# =============================================================================
# 06_ncycling_guilds.R
# MICP Microbiome Analysis Pipeline — N-Cycling Functional Guild Analysis
# =============================================================================
# PURPOSE:
#   Assigns OTUs to nitrogen-cycling functional guilds based on genus-level
#   taxonomy. Computes guild-level summed relative abundance per sample,
#   tests for treatment effects (Kruskal-Wallis + Dunn), and generates
#   faceted boxplots.
#
# INPUT:  data/processed/ps.rds  (full, non-rarefied object)
# OUTPUT: outputs/figures/fig_ncycling_guilds_by_treatment.{tiff,png}
#         outputs/tables/ncycling_guild_summary.csv
#         outputs/tables/ncycling_guild_stats.csv
#
# DEPENDENCIES: tidyverse, phyloseq, rstatix
#
# NOTES ON GUILD ASSIGNMENT:
#   Guild definitions are maintained in N_GUILDS (00_config.R). Genera with
#   documented roles in multiple N-transformation pathways (e.g., Paenibacillus
#   as Ureolytic, Denitrifier, and Diazotroph) contribute their relative
#   abundance independently to each applicable guild. Guild-level relative
#   abundances are therefore NOT mutually exclusive and do NOT sum to unity
#   across guilds. This is by design and should be stated explicitly in any
#   manuscript methods section.
#
#   Relevance to MICP sustainability:
#   Ureolytic guild tracks taxa that could contribute to or compete with
#   S. pasteurii-driven ureolysis. AOB/NOB track nitrification capacity,
#   relevant to NH4+ fate following MICP treatment. Denitrifiers and
#   diazotrophs reflect the broader N-cycling ecosystem services potentially
#   disrupted by the elevated NH4+ and altered pH from MICP treatment.
# =============================================================================

source("R/00_config.R")

library(tidyverse)
library(phyloseq)
library(rstatix)


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
# 1. LOAD & PREPARE
# -----------------------------------------------------------------------------
cat("[06] Loading phyloseq object...\n")
ps <- readRDS("data/processed/ps.rds")

# Build relative abundance object (non-rarefied; relative abundance normalizes
# for library size differences in this context)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))


# -----------------------------------------------------------------------------
# 2. BUILD OTU → GUILD MAP
# -----------------------------------------------------------------------------
tax_df      <- as.data.frame(tax_table(ps))
tax_df$OTU  <- rownames(tax_df)

guild_map <- map_dfr(names(N_GUILDS), function(guild_name) {
  matched_otus <- tax_df$OTU[tax_df$GENUS %in% N_GUILDS[[guild_name]]]
  if (length(matched_otus) == 0) return(NULL)
  data.frame(OTU = matched_otus, guild = guild_name, stringsAsFactors = FALSE)
})

# Diagnostic: report OTUs assigned to multiple guilds
multi_guild <- guild_map %>% count(OTU) %>% filter(n > 1)
cat(sprintf("[06] OTUs assigned to multiple guilds (expected): %d\n",
            nrow(multi_guild)))

# Report coverage per guild
guild_coverage <- guild_map %>%
  count(guild, name = "n_OTUs") %>%
  arrange(desc(n_OTUs))
cat("\nGuild coverage (OTUs matched):\n")
print(guild_coverage)

# Warn if any defined guild matched nothing
missing_guilds <- setdiff(names(N_GUILDS), guild_coverage$guild)
if (length(missing_guilds) > 0) {
  warning(sprintf(
    "[06] No OTUs matched for guild(s): %s\n  Check genus names in N_GUILDS (00_config.R) against taxonomy in your dataset.",
    paste(missing_guilds, collapse = ", ")
  ))
}


# -----------------------------------------------------------------------------
# 3. COMPUTE GUILD-LEVEL RELATIVE ABUNDANCE PER SAMPLE
# -----------------------------------------------------------------------------
otu_rel <- as.data.frame(otu_table(ps_rel))
if (!taxa_are_rows(ps_rel)) otu_rel <- as.data.frame(t(otu_rel))
otu_rel$OTU <- rownames(otu_rel)

otu_long <- otu_rel %>%
  pivot_longer(-OTU, names_to = "sample", values_to = "rel_abund")

# One-to-many join: multi-guild OTUs contribute independently to each guild
guild_summary <- otu_long %>%
  inner_join(guild_map, by = "OTU") %>%
  group_by(sample, guild) %>%
  summarise(total_rel_abund = sum(rel_abund), .groups = "drop")

# Attach treatment metadata
meta_df        <- as.data.frame(sample_data(ps))
meta_df$sample <- rownames(meta_df)

guild_summary <- guild_summary %>%
  left_join(meta_df[, c("sample", COL_NAMES$treatment)], by = "sample") %>%
  rename(Treatment = all_of(COL_NAMES$treatment)) %>%
  mutate(Treatment = factor(Treatment, levels = TREATMENT_MAP$display_labels),
         guild     = factor(guild, levels = names(N_GUILDS)))


# -----------------------------------------------------------------------------
# 4. STATISTICS — Kruskal-Wallis + Dunn pairwise per guild
# -----------------------------------------------------------------------------
cat("\n=== Guild Analysis: Kruskal-Wallis Tests ===\n")
guild_kw <- guild_summary %>%
  group_by(guild) %>%
  kruskal_test(total_rel_abund ~ Treatment) %>%
  adjust_pvalue(method = GUILD_POSTHOC_METHOD)
print(guild_kw)

cat(sprintf("\n=== Pairwise Dunn Tests (%s correction) ===\n",
            GUILD_POSTHOC_METHOD))
guild_dunn <- guild_summary %>%
  group_by(guild) %>%
  dunn_test(total_rel_abund ~ Treatment, p.adjust.method = GUILD_POSTHOC_METHOD)
print(guild_dunn)

# Save stats
write.csv(guild_kw,
          file.path(OUT_DIR_TABLES, "ncycling_guild_stats_kw.csv"),
          row.names = FALSE)
write.csv(guild_dunn,
          file.path(OUT_DIR_TABLES, "ncycling_guild_stats_dunn.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 5. PLOT
# -----------------------------------------------------------------------------
# Dynamic ncol: use 3 columns for ≤6 guilds, else auto
n_guilds_plot <- length(unique(guild_summary$guild))
ncol_guilds   <- min(3, n_guilds_plot)

p_guild <- ggplot(guild_summary,
                  aes(x = Treatment, y = total_rel_abund, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  facet_wrap(~ guild, scales = "free_y", ncol = ncol_guilds) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title   = "N-cycling guild representation by MICP treatment",
    x       = "Treatment",
    y       = "Summed relative abundance",
    caption = "Multi-guild membership permitted; abundances are not mutually exclusive across guilds."
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold"),
    axis.text.x     = element_text(angle = 35, hjust = 1)
  )

print(p_guild)
save_figure(p_guild, "fig_ncycling_guilds_by_treatment", width = 10, height = 6)


# -----------------------------------------------------------------------------
# 6. EXPORT GUILD SUMMARY TABLE
# -----------------------------------------------------------------------------
write.csv(guild_summary,
          file.path(OUT_DIR_TABLES, "ncycling_guild_summary.csv"),
          row.names = FALSE)
cat(sprintf("[06] Guild summary table saved to %s/ncycling_guild_summary.csv\n",
            OUT_DIR_TABLES))

cat("[06] N-cycling guild analysis complete.\n")
