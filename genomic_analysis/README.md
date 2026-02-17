# Genomic Analysis Pipeline

This directory contains scripts for whole-genome sequencing (WGS) analysis of evolved clones: alignment, SNP calling, coverage-based structural variant detection, and barcode QC.

## Pipeline overview

### 1. (Optional) WGS preprocessing — `wgs_preprocessing/`

**Script:** `merge_wgs.py`

Some clones were sequenced twice to achieve higher coverage. This script concatenates the FASTQ.gz files from two sequencing runs for the same clone using `cat` (preserving gzip format).

Required inputs:
- `../filename_info/resequenced_wgs.csv` — maps clone indices to file keys from each sequencing run
- Raw FASTQ.gz files from two sequencing runs

Output: merged FASTQ.gz files in `merged_files/`

---

### 2. Alignment and variant calling — `snp_calling/`

**Script:** `alignment.sh`

Submitted as a SLURM array job; each task processes one sample.

1. **Align** paired-end reads to the S288C R64-3-1 reference genome with `bwa mem`.
2. **Add read groups** with Picard `AddOrReplaceReadGroups`.
3. **Mark duplicates** with Picard `MarkDuplicates`.
4. **Call variants** per-sample in GVCF mode with GATK `HaplotypeCaller`.

Outputs: per-sample `.sorted.rg.md.bam` and `.g.vcf.gz` files.

---

### 3. Joint genotyping — `snp_calling/`

**Script:** `merge_vcf.sh`

Run after `alignment.sh` has finished for all samples.

1. **Consolidate GVCFs** across all samples with GATK `GenomicsDBImport`, using `intervals.list` to specify chromosomes and `sample_map.txt` to specify samples.
2. **Joint genotype** with GATK `GenotypeGVCFs`.

**`sample_map.txt` format** (not included; contains cluster-specific absolute paths):
A tab-separated file with two columns: sample name and absolute path to the sample's `.g.vcf.gz` file. One sample per line.

Output: `wgs_all_samples.vcf.gz`

---

### 4. SNP filtering — `snp_calling/`

**Script:** `basic_filter.py`

Applies hard quality filters to the joint-genotyped VCF and outputs candidate SNPs as a CSV.

Filter thresholds (stricter than GATK defaults):
| Annotation | Direction | Threshold |
|---|---|---|
| QD (quality by depth) | > | 4 |
| FS (Fisher strand bias) | < | 40 |
| SOR (strand odds ratio) | < | 2 |
| MQ (mapping quality) | > | 50 |
| MQRankSum | > | -3 |
| ReadPosRankSum | > | -5 |

Additionally filters out variants with insufficient alternate allele depth (< 2 reads) or called identically in more than 3 samples (likely ancestral polymorphisms).

Outputs:
- `../bamsnap/snp_coords.csv` — per-SNP table with sample, genomic location, gene name, evolutionary condition
- `../bamsnap/dup_snp_coords.csv` — SNPs called in more than one sample (for manual review)
- `filtered_snps.vcf.gz` — VCF with passing SNPs only

---

### 5. SNP validation — `snp_calling/`

**Script:** `gen_snaps.py`

Generates bamsnap images for visual inspection of each candidate SNP. Each image shows read alignments, per-base identity, and coverage at the variant site, with an ancestral control sample included for comparison.

Usage: `python gen_snaps.py <snp_coords_csv> <output_dir> [num_snaps]`

---

### 6. Coverage analysis — `coverage_analysis/`

Used to detect large-scale copy number variants (CNVs) and aneuploidies.

**Step 6a — Generate coverage** (`gen_coverage.py`): Runs `samtools depth -a` on each BAM to produce a per-position coverage file. Submitted per-sample (indexed by line in `wgs_filenames.txt`).

Usage: `python gen_coverage.py <genome_index>`

**Step 6b — Coverage by gene** (`cov_by_gene.py`): Reads coverage files and computes mean coverage per annotated gene using the S288C R64-3-1 GFF. Saves output as `cov_by_gene_data.npz`.

**Step 6c — Regional coverage** (`regional_cov.py`): Computes windowed mean coverage in a user-specified genomic region (default: chrIII:200000-250000, 100 bp bins). Edit `chrom_id`, `start_pos`, `end_pos`, and `bin_size` at the top of the script as needed. Saves output as `regional_cov.npz`.

---

### 7. QC — Barcode recovery from WGS — `barcode_from_wgs/`

A quality control step to verify that WGS reads for each clone span and capture the expected barcode sequence.

**Step 7a — Align to barcode locus** (`align_bc.sh`): SLURM array job that aligns WGS reads to a reference FASTA containing only the barcode insertion locus (`YBR209W_locus`).

**Step 7b — Extract barcode** (`extract_read_info.py`): Reads the resulting BAM files and recovers the 26 bp barcode sequence by majority vote across all reads covering that position, including soft-clipped bases. Outputs `retrieved_barcodes.csv` for comparison against the expected barcode assignments.
