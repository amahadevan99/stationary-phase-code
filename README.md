# Stationary Phase Evolution — Processing Scripts

This repository contains the data processing scripts used in our study of adaptive evolution in *Saccharomyces cerevisiae* propagated in non-fermentable (stationary phase) carbon sources. The paper describes bulk fitness assays using random-barcode lineage tracking and whole-genome sequencing of evolved clones.

## Repository layout

```
stationary_phase_repo/
├── figures_data_code/          # Data and scripts needed to generate all paper figures
│   ├── R_scripts/              # R scripts for figure generation
│   ├── python_scripts/         # Python scripts for data processing and figures
│   ├── R_data/                 # Final fitness data with genotypes attached (input to R scripts)
│   ├── python_data/            # Processed data files (input to Python scripts)
│   ├── fitseq2_input/          # FitSeq2 input count files
│   ├── fitseq2_output/         # FitSeq2 results
│   ├── fitmut_input/           # FitMut input count files (filtered, unfiltered, params)
│   ├── fitmut_output_filtered/ # FitMut results
│   └── T0_counts_aggregated/   # Aggregated T0 counts
├── barcode_analysis/
│   ├── barcode_counting/       # Demultiplex reads, extract & cluster barcodes; map clones to fitness
│   ├── fitness_remeasurement/  # Independent fitness remeasurement pipeline for selected clones
│   └── fitness_inference/      # Reformat barcode counts; run FitSeq2 / FitMut2
└── genomic_analysis/
    ├── wgs_preprocessing/      # Concatenate FASTQ from multiple sequencing runs
    ├── snp_calling/            # Alignment, joint genotyping, hard filtering, visual validation
    ├── coverage_analysis/      # Per-position and per-gene coverage; structural variant detection
    └── barcode_from_wgs/       # QC: recover barcode from WGS reads
```

See `figures_data_code/data_files_README.md` for detailed descriptions of the data files.

See `barcode_analysis/README.md` and `genomic_analysis/README.md` for step-by-step pipeline descriptions.

## Raw data

Raw genomic sequencing data are deposited at [external archive — DOI TBD].

## Software dependencies

| Tool | Purpose |
|---|---|
| [bartender](https://github.com/LaoZZZZZ/bartender-1.1) | Barcode extraction and clustering |
| [FitSeq2](https://github.com/FangfeiLi05/FitSeq2) | Fitness inference from fitness measurement trajectories with adaptive mutants |
| [FitMut2](https://github.com/FangfeiLi05/FitMut2) | Fitness inference from experimental evolutions |
| [BWA](https://github.com/lh3/bwa) | Short-read alignment to reference genome |
| [Samtools](http://www.htslib.org/) | BAM manipulation and depth calculation |
| [Picard](https://broadinstitute.github.io/picard/) | Read group annotation and duplicate marking |
| [GATK](https://gatk.broadinstitute.org/) | Variant calling and joint genotyping |
| [bamsnap](https://github.com/parklab/bamsnap) | Visual SNP validation from BAM files |

Python package dependencies include: `numpy`, `pandas`, `cyvcf2`, `pysam`, `biopython`, `scipy`, `python-Levenshtein`.

## Cluster notes

Shell scripts (`.sh`) are written as SLURM array jobs for the [Sherlock HPC cluster](https://www.sherlock.stanford.edu/) at Stanford. Absolute paths to reference genomes, software, and scratch storage are cluster-specific and will need to be updated for other environments.

Thanks to Claude Code for documenting and organizing this repository.
