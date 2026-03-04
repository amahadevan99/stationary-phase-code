# Data Files

## R_data/

Final fitness and performance data, after correcting for index-hopping, with all genotype data attached.

- **adaptive_haploid_lineages.csv** — Only adaptive haploid lineages are included.
- **adaptive_haploid_no_SMF2.csv** — Only adaptive haploid lineages without an SMF2 mutation are included.
- **adaptive_haploid_SMF2.csv** — Only adaptive haploid SMF2 mutations are included.
- **adaptive_lineages.csv** — Only adaptive lineages are included, except for 1/5 day evolved which are all included.

## python_data/

- **adaptive_fitness_home_condition.csv** — Only adaptive lineages are included. No 1/5 day lineages included.
- **lane1_raw_counts.csv** and **lane2_raw_counts.csv** — Raw counts from the two sequencing lanes for the fitness measurement experiments.
- **lane1_chimera_corrected.csv** and **lane2_chimera_corrected.csv** — Counts after correcting for index hopping.
- **fitseq2_fitness_all_replicates.csv** — FitSeq2 outputs (from `fitseq2_output/`) for the 27 fitness measurement experiments arrayed in a single CSV file.
- **chimeric_predictions_empty_combos.csv** — Real counts and predicted counts for primer combos for which no sample was prepared.
- **fitness_data_comparison.csv** — Data for clones measured previously, containing this study's data alongside Li et al. and Kinsler et al. data.
- **evo_transfer_record.csv** — Fraction of clones that were diploid at each transfer during the evolution.
- **all_samples.csv** — All fitness measurement BC sequencing samples, listing the ROS and FOS primers used, the sequencing lane, and whether the same timepoint was sequenced on both lanes.

## fitseq2_input/

Index-hopping corrected and ancestor-aggregated BC counts used as input to FitSeq2. Contains 27 CSV files, one per experiment (e.g., `M3_1_day_r1.csv`).

## fitseq2_output/

FitSeq2 output from the input above. Contains 3 files per experiment (81 total): `*_FitSeq2_Result.csv`, `*_Mean_fitness.csv`, and `*_Read_Number_Estimated.csv`.

## fitmut_input/

Input data for FitMut. Contains three subdirectories:

- **filtered/** — Filtered barcode count matrices (9 CSV files, one per experiment), used as input to FitMut after removing uneven barcode trajectories.
- **unfiltered/** — Unfiltered barcode count matrices (9 CSV files).
- **params0/** — FitMut parameter files including `*_timepoints.csv` (generation numbers and population sizes for each timepoint).

## fitmut_output_filtered/

FitMut output files. Contains `*_Cell_Number.csv` and `*_MutSeq_Result.csv` pairs for each experiment.

## T0_counts_aggregated/

Aggregated T0 (timepoint zero) counts used to prepend to trajectory data: `FOS3_ROS1_agg.csv`, `FOS5_ROS9_agg.csv`, `FOS9_ROS7_agg.csv`.
