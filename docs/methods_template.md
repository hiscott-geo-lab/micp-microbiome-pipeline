# Methods Template — Microbiome Data Analysis

*Copy-paste into manuscript methods section and fill in bracketed values from your run. All statistical values should be verified against pipeline console output and exported tables.*

---

## Bioinformatics & Community Analysis

Raw Nanopore reads were processed using NanoClust and NanoScript, Nextflow-based pipelines for long-read microbiome analysis. Taxonomic placement was performed using the T-BAS platform (v2.4) with the LIFE1 reference phylogeny, yielding [X] OTUs across [N] samples. Samples from [EXCLUDED_SPECIES/CONDITION] were removed prior to analysis [state reason].

Downstream community analyses were conducted in R (v[X.X.X]) using the phyloseq package (v[X.X.X]; McMurdie & Holmes 2013). MICP treatment labels are defined as follows: U250 (equimolar urea and CaCl₂ at 0.25 M), U100+P (0.1 M urea-CaCl₂ with nutrient amendment), and U50+P (0.05 M urea-CaCl₂ with nutrient amendment).

### Relative Abundance

Sequences were agglomerated at the genus level and converted to relative abundance. The [TOP_N_GENERA] most abundant genera across all samples are shown in stacked bar plots averaged by treatment and faceted by plant species.

### Alpha Diversity

Alpha diversity was assessed using Shannon entropy, Simpson's index, observed OTU richness, and Faith's phylogenetic diversity (Faith's PD; Faith 1992), computed using the picante package (v[X.X.X]). Note that singletons were removed during the T-BAS pipeline; observed OTU richness values are therefore conservative and should be interpreted in relative rather than absolute terms. Treatment effects on alpha diversity were assessed using Kruskal-Wallis tests, with pairwise post-hoc comparisons performed using Dunn's test with [BONFERRONI/BH] correction.

### Beta Diversity

Beta diversity was assessed using Bray-Curtis dissimilarity and unweighted and weighted UniFrac distances (Lozupone & Knight 2005). Prior to ordination, samples were rarefied to [RAREFY_DEPTH] reads (the [minimum observed/user-specified] sequencing depth) using a fixed random seed for reproducibility. Community composition was visualized by Principal Coordinates Analysis (PCoA) and Non-metric Multidimensional Scaling (NMDS; stress = [X.XXX]).

The effect of MICP treatment and plant species on community composition was tested using PERMANOVA (Anderson 2001) via the adonis2 function in vegan (v[X.X.X]; Oksanen et al. [YEAR]), with [PERMANOVA_PERMS] permutations and sequential ("by = terms") partitioning. Pairwise treatment contrasts were assessed using PERMANOVA on distance subsets, with Bonferroni correction for multiple comparisons. Homogeneity of multivariate dispersions was tested using PERMDISP (betadisper; Anderson 2006) as a validity check for PERMANOVA interpretation.

### N-Cycling Functional Guilds

OTUs were assigned to nitrogen-cycling functional guilds (ureolytic, ammonia-oxidizing bacteria [AOB], nitrite-oxidizing bacteria [NOB], denitrifiers, diazotrophs) based on genus-level taxonomy and published metabolic annotations. Multi-guild membership was permitted for genera with documented roles in multiple N-transformation pathways; guild-level relative abundances are therefore non-exclusive and do not sum to unity across guilds. Guild-level relative abundances were compared across treatments using Kruskal-Wallis tests with Dunn post-hoc tests ([BH/BONFERRONI] correction).

---

## Key References

- Anderson, M.J. (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecology*, 26, 32–46.
- Anderson, M.J. (2006). Distance-based tests for homogeneity of multivariate dispersions. *Biometrics*, 62, 245–253.
- Faith, D.P. (1992). Conservation evaluation and phylogenetic diversity. *Biological Conservation*, 61, 1–10.
- Lozupone, C. & Knight, R. (2005). UniFrac: a new phylogenetic method for comparing microbial communities. *Applied and Environmental Microbiology*, 71, 8228–8235.
- McMurdie, P.J. & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE*, 8, e61217.
- Oksanen, J. et al. ([YEAR]). vegan: Community Ecology Package. R package version [X.X.X].
