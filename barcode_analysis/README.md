# Barcode Analysis Pipeline

This directory contains scripts for tracking lineage fitness through barcode sequencing (BarSeq) and for associating WGS-sequenced clones with their inferred fitness values.

## Pipeline overview

### 1. Barcode counting — `barcode_counting/`

**Script:** `barcode_counts.py`

Processes paired-end barcode sequencing data to produce per-lineage read-count trajectories.

1. **Demultiplex** raw FASTQ reads by first-step PCR primer identity (forward/reverse multitags).
2. **Extract barcodes** from demultiplexed reads using a flanking sequence anchor (`TATCGGTACC`); optionally extract UMIs.
3. **Cluster barcodes** per sample/timepoint using `bartender_single_com` (edit distance ≤ 3).
4. **Combine timepoints** across sequencing dates using `bartender_combiner_com` to produce one cluster file per evolution experiment.

Required inputs:
- `{date}_barcode_files.txt` — list of raw FASTQ filenames for that sequencing date
- `barcode_info/{date}/multitags_list.csv` — primer tag → name mapping
- `barcode_info/{date}/filename_key.csv` — filename → N,S primer name mapping
- `barcode_info/{date}/bc_indices.csv` — timepoint label → primer combination mapping

Outputs:
- `combined_bc1_find_loc/{exp_id}_cluster.csv` — bartender combined cluster output (read counts per timepoint per lineage)

---

### 2. Fitness inference — `fitness_inference/`

**Scripts:** `reformat_bc_counts.py` → `fitseq_inference.py` (primary) or `fitmut_inference.py`

**Step 2a — Reformat** (`reformat_bc_counts.py`): Reads bartender cluster output, filters lineages with low inverse participation ratio (IPR < 1.3 by default), and writes a plain count matrix for use by inference tools.

**Step 2b — Infer fitness** (`fitseq_inference.py`): Runs FitSeq2 on the reformatted count matrix. FitSeq2 is the primary inference method used in the paper. Accepts an experiment index as a command-line argument, which selects from the hardcoded experiment list.

`fitmut_inference.py` runs FitMut2 as an alternative/validation method.

Required inputs for FitSeq2:
- `fitseq_input/{exp_name}_traj.csv` — barcode count trajectories
- `fitseq_input/{exp_name}_timepoints.csv` — timepoints in generations

---

### 3. Clone identification — `barcode_counting/`

**Script:** `clone_fitnesses.py`

Associates WGS-sequenced evolved clones with their barcode lineage and inferred fitness.

1. **Extract barcodes** from clone barcode sequencing FASTQ, demultiplexed by plate-specific primer pairs.
2. **Match** each clone's majority barcode sequence against the bartender cluster centers.
3. **Assign fitness** from FitMut2/FitSeq2 output for the matching barcode cluster.

Required inputs:
- Clone barcode sequencing FASTQ (paired-end)
- `clone_info/forward_primers.csv`, `clone_info/reverse_primers.csv` — primer tag sequences
- `clone_info/{date}_clone_info.csv` — clone metadata
- Bartender cluster output and fitness inference results

Outputs:
- `{date}/barcodes.csv` — per-clone majority barcode(s) and read counts
- `{date}/labeled_barcodes.csv` — barcodes matched to cluster with inferred fitness
- `{date}/unique_barcodes.csv` — deduplicated clone barcode table

---

### 4. Fitness remeasurement — `fitness_remeasurement/`

**Scripts:** `extract_bcs.py` → `count_bcs.py`

A separate, simpler pipeline for remeasured clones that uses a known reference barcode list rather than de novo clustering.

**Step 4a — Extract** (`extract_bcs.py`): Demultiplexes reads by FOS/ROS primer pairs (Hamming distance ≤ 2) and extracts barcode sequences using a flanking sequence anchor.

Usage: `python extract_bcs.py <date>`

**Step 4b — Count** (`count_bcs.py`): For each demultiplexed primer-pair file, counts occurrences of each barcode in the reference list (with ±1 bp length tolerance).

Usage: `python count_bcs.py <date>`

Required inputs:
- `{date}_barcode_files.txt` — list of raw FASTQ filenames
- `barcode_info/primer_indices.csv` — FOS/ROS primer tag sequences
- `barcode_info/barcode_list.csv` — reference barcode sequences to count

Output:
- `bc_counts_{date}.csv` — count matrix (rows = primer pairs/clones, columns = barcodes)
