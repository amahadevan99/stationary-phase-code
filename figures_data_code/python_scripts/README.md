# Python Scripts for Figure Generation

All scripts must be run from the `python_scripts/` directory. Output PDFs are saved to `python_scripts/figures/`.

## Dependencies

```
pandas numpy matplotlib seaborn scikit-learn scipy
```

## Scripts

### fig_2a_entropy.py — Figure 2a
Plots Shannon diversity of barcode frequencies over generations for each Gly/Eth evolution condition. Lines are colored by transfer interval.

**Input:**
- `../fitmut_input/filtered/*_counts.csv`
- `../fitmut_input/params0/*_timepoints.csv`

### fig_2b_fitness_effect.py — Figure 2b
Plots fitness effect per cycle in the evolution condition vs. days between transfers, for Glucose and Gly/Eth media. Shows jittered scatter with mean bars and linear regression.

**Input:**
- `../python_data/adaptive_fitness_home_condition.csv`

### fig_S1_trajectories.py — Supplemental Figure S1
Plots lineage frequency trajectories (top 2,000 by total count) from fitmut output, colored by fitness. One subplot per evolution condition/replicate, with star markers indicating specific timepoints.

**Input:**
- `../fitmut_output_filtered/*_Cell_Number.csv` (one per condition/replicate)
- `../fitmut_output_filtered/*_MutSeq_Result.csv` (matching fitness results)

### fig_S2_top100_trajectories.py — Supplemental Figure S2
Plots measured frequency trajectories for the top 100 fittest lineages plus 100 random others, using FitSeq2 output. Generates separate Glucose (M3) and Gly/Eth (M05) panel PDFs. Prepends T0 counts and annotates total read counts per timepoint.

**Input:**
- `../fitseq2_input/*.csv` (trajectory files, e.g. `M3_1_day_r1.csv`)
- `../T0_counts_aggregated/FOS3_ROS1_agg.csv`
- `../T0_counts_aggregated/FOS5_ROS9_agg.csv`
- `../T0_counts_aggregated/FOS9_ROS7_agg.csv`
- `../fitseq2_output/*_FitSeq2_Result.csv`

### fig_S3_replicate_correlation.py — Supplemental Figure S3
Replicate-vs-replicate fitness correlation plots (R1 vs R2, R1 vs R3, R2 vs R3) for each evolution condition. Generates separate PDFs for Glucose and Gly/Eth. Filters lineages with error > 5.

**Input:**
- `../python_data/fitseq2_fitness_all_replicates.csv`

### fig_S4_method_comparison.py — Supplemental Figure S4
Pairwise comparison of fitness measurements across three methods (Li et al. 2019, Kinsler et al. 2024, this study) at 1-, 2-, 3-, and 5-day transfer intervals. Scatter plots colored by ancestor background with R² statistics.

**Input:**
- `../python_data/fitness_data_comparison.csv`

### fig_S13_growth_curve.py — Supplemental Figure S13
Growth curve scatter plot showing CFU/mL (from plating) over elapsed time in hours, with a vertical line at 48 hours.

**Input:**
- `../python_data/growth_curve_data.csv`

### fig_S14_ploidy.py — Supplemental Figure S14
Plots percent diploid over generations for each Gly/Eth evolution condition (2-, 4-, 6-, 8-, 10-day transfers), with separate lines per replicate.

**Input:**
- `../python_data/evo_transfer_record.csv`

### fig_S15_chimeric_predictions.py — Supplemental Figure S15
Log-log scatter plots of predicted vs. observed chimeric read counts in empty ROS/FOS combinations. Shows regression line and x=y reference.

**Input:**
- `../python_data/chimeric_predictions_empty_combos.csv`

### fig_S16_PCR_replicates.py — Supplemental Figure S16
PCR replicate correlation plots comparing raw counts and chimera-corrected counts across sequencing lanes. Pairs are sorted by media type and transfer interval.

**Input:**
- `../python_data/lane1_raw_counts.csv`
- `../python_data/lane2_raw_counts.csv`
- `../python_data/lane1_chimera_corrected.csv`
- `../python_data/lane2_chimera_corrected.csv`
- `../python_data/all_samples.csv`
