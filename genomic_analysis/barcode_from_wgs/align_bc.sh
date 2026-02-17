#!/bin/bash
# Aligns WGS reads to the barcode insertion locus reference for QC verification.
#
# Usage:
#   Submitted as a SLURM array job on the Sherlock HPC cluster.
#   Each array task processes one sample, identified by line number in wgs_filenames.txt.
#   Adjust --array range to match the number of samples.
#
# Inputs:
#   - wgs_filenames.txt                     — one sample filename prefix per line
#   - {dir}{bn}_R1_001.fastq.gz / _R2_001.fastq.gz  — paired-end WGS reads
#   - barcode_region/barcode_region.fasta   — reference FASTA for the barcode insertion locus
#
# Outputs:
#   - bc_bam/{sample}.sorted.bam      — sorted BAM aligned to barcode locus
#   - bc_bam/{sample}.sorted.bam.bai  — BAM index
#
# Notes:
#   This is a QC step to verify that WGS reads capture the expected barcode sequences.
#   Use extract_read_info.py to retrieve the barcode sequence from the resulting BAM.
#
# Dependencies:
#   BWA, samtools; loaded via `ml restore bioinformatics`
#SBATCH --job-name=alignment
#SBATCH --time=48:00:00
#SBATCH -p normal
#SBATCH --cpus-per-task=4
#SBATCH --array=11-576
#SBATCH --output=slurm_files/alignment_%a.txt

ml restore bioinformatics

bn=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/groups/dsfisher/adityam/sherlock_lab/cerevisiae/genomes/wgs_filenames.txt)

ref_genome=$"/home/groups/dsfisher/adityam/sherlock_lab/cerevisiae/genomes/align_to_barcode/barcode_region/barcode_region.fasta"
dir=$"/scratch/users/adityam2/wgs_data/"
out_dir=$"/scratch/users/adityam2/wgs_data/bc_bam/"

### generates bam files; puts output files in the same directory regardless of which sequencing run they were from
bwa mem -t 4 "$ref_genome" "$dir$bn"_R1_001.fastq.gz "$dir$bn"_R2_001.fastq.gz | samtools view -S -b - | samtools sort > "$out_dir${bn##*/}".sorted.bam
samtools index "$out_dir${bn##*/}".sorted.bam
