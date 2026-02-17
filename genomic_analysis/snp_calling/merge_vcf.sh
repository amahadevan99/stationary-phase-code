#!/bin/bash
# Merges per-sample GVCFs and performs joint genotyping across all samples.
#
# Usage:
#   Submitted as a SLURM job on the Sherlock HPC cluster.
#   Run after alignment.sh has produced per-sample GVCFs for all samples.
#
# Inputs:
#   - intervals.list      — chromosome list for GATK (S288C R64-3-1 chromosomes)
#   - sample_map.txt      — tab-separated file mapping sample names to GVCF paths
#                           (format: <sample_name>\t<path_to_gvcf>)
#   - S288C R64-3-1 reference genome (.fasta)
#
# Outputs:
#   - {out_dir}genomicsdb/  — GenomicsDB workspace (intermediate)
#   - wgs_all_samples.vcf.gz  — joint-genotyped VCF across all samples
#
# Dependencies:
#   GATK; loaded via `ml restore bioinformatics`
#SBATCH --job-name=merge_vcf
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH -p normal
#SBATCH --output=slurm_files/merge_vcf.txt

ml restore bioinformatics

ref_genome=$"/home/groups/dsfisher/adityam/sherlock_lab/cerevisiae/S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fasta"
out_dir=$"/scratch/users/adityam2/wgs_data/"

gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport -L intervals.list --sample-name-map sample_map.txt --genomicsdb-workspace-path "$out_dir"genomicsdb

gatk --java-options "-Xmx4g" GenotypeGVCFs -R "$ref_genome" -V gendb://"$out_dir"/genomicsdb -O wgs_all_samples.vcf.gz

# bcftools view -S 10X_cov_samples.txt wgs_all_samples.vcf.gz > wgs_all_samples_10X_coverage.vcf

# gzip wgs_all_samples_10X_coverage.vcf
