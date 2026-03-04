# R Scripts for Figure Generation

All scripts must be run from the `R_scripts/` directory. Output PDFs are saved to `R_scripts/figures/`.

## Dependencies

```
ggplot2 dplyr broom grid stringr tidyr
```

## Scripts

### figure_3a.R — Figure 3a
Scatter plot of earliest stationary phase vs. glucose stationary phase fitness change from ancestor, colored by ancestor and shaped by ploidy. Linear regression per ancestor group with R² statistics.

**Input:** `../R_data/adaptive_lineages.csv`

### figure_3b.R — Figure 3b
Faceted scatter plot (3x3) of earliest stationary phase vs. glucose stationary phase fitness change from ancestor, per evolution condition. Points colored by ancestor, with linear regression and R²/p-value per panel.

**Input:** `../R_data/adaptive_haploid_lineages.csv`

### figure_4a.R — Figure 4a
Scatter plot of late vs. earliest stationary phase fitness change from ancestor, colored by ancestor and shaped by ploidy. Linear regression per ancestor group with R² statistics.

**Input:** `../R_data/adaptive_lineages.csv`

### figure_4b.R — Figure 4b
Faceted scatter plot (3x3) of late vs. earliest stationary phase fitness change from ancestor, per evolution condition. Points colored by ancestor with linear regression lines (Gly/Eth 8 day excluded from trendlines).

**Input:** `../R_data/adaptive_haploid_lineages.csv`

### figure_5.R — Figure 5
Scatter plot of late vs. earliest stationary phase fitness change for SMF2-containing lineages, colored by ancestor and shaped by ploidy. Linear regression per ancestor with R² statistics.

**Input:** `../R_data/adaptive_haploid_SMF2.csv`

### figure_S5.R — Supplemental Figure S5
Jitter plot of fitness change per cycle in the Gly/Eth 6-day condition, grouped by genotype (chr11 duplication + SMF2, chr11 duplication only, SMF2 only).

**Input:** `../R_data/adaptive_lineages.csv`

### figure_S6.R — Supplemental Figure S6
Jitter plot of fitness change per cycle across all assay environments, grouped by FZF1 diploid, diploid, and haploid genotypes. Points sized by group.

**Input:** `../R_data/adaptive_lineages.csv`

### figure_S7.R — Supplemental Figure S7
Faceted scatter plot (3x3) of earliest stationary phase vs. glucose stationary phase fitness change from ancestor, per evolution condition. Same as Figure 3b but includes all ploidy states (haploid, diploid, undetermined).

**Input:** `../R_data/adaptive_lineages.csv`

### figure_S8.R — Supplemental Figure S8
Scatter plot of late vs. earliest stationary phase fitness change for the Gly/Eth 6-day condition only, colored by SMF2 mutation status. Linear regression per SMF2 group.

**Input:** `../R_data/adaptive_haploid_lineages.csv`

### figure_S9a.R — Supplemental Figure S9a
Scatter plot of latest vs. early stationary phase fitness change from ancestor, colored by ancestor and shaped by ploidy. Trendlines exclude Levy et al. 2015 ancestor.

**Input:** `../R_data/adaptive_lineages.csv`

### figure_S9b.R — Supplemental Figure S9b
Faceted scatter plot (3x3) of latest vs. early stationary phase fitness change from ancestor, per evolution condition. Points colored by ancestor with linear regression.

**Input:** `../R_data/adaptive_haploid_lineages.csv`

### figure_S10a.R — Supplemental Figure S10a
Scatter plot of late vs. earliest stationary phase fitness change from ancestor (SMF2 lineages excluded), colored by ancestor and shaped by ploidy.

**Input:** `../R_data/adaptive_haploid_no_SMF2.csv`

### figure_S10b.R — Supplemental Figure S10b
Scatter plot of latest vs. early stationary phase fitness change from ancestor (SMF2 lineages excluded), colored by ancestor and shaped by ploidy.

**Input:** `../R_data/adaptive_haploid_no_SMF2.csv`

### figure_S11.R — Supplemental Figure S11
Scatter plot of latest vs. early stationary phase fitness change for SMF2-containing lineages, colored by ancestor and shaped by ploidy.

**Input:** `../R_data/adaptive_haploid_SMF2.csv`

### figure_S12.R — Supplemental Figure S12
Violin + boxplot + jitter of fitness change per cycle for chr11 duplication lineages (excluding FZF1) across assay environments.

**Input:** `../R_data/adaptive_haploid_lineages.csv`
