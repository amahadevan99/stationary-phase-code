#!/bin/bash
# Aligns WGS reads to the S288C reference genome, marks duplicates, and calls variants (GVCF mode).
#
# Usage:
#   Submitted as a SLURM array job on the Sherlock HPC cluster.
#   Each array task processes one sample, identified by line number in wgs_filenames.txt.
#   Adjust --array range to match the number of samples.
#
# Inputs:
#   - wgs_filenames.txt             — one sample filename prefix per line (without _R1/_R2 suffix)
#   - {dir}{bn}_R1_001.fastq.gz / _R2_001.fastq.gz  — paired-end WGS reads
#   - S288C R64-3-1 reference genome (.fsa for BWA, .fasta for GATK)
#
# Outputs:
#   - bam/{sample}.sorted.bam           — sorted BAM
#   - bam/{sample}.sorted.rg.bam        — BAM with read groups added
#   - bam/{sample}.sorted.rg.md.bam     — BAM with duplicates marked (used for downstream steps)
#   - vcf/{sample}.g.vcf.gz             — per-sample GVCF for joint genotyping
#
# Dependencies:
#   BWA, samtools, Picard (picard.jar), GATK; loaded via `ml restore bioinformatics`
#SBATCH --job-name=alignment
#SBATCH --time=48:00:00
#SBATCH -p normal
#SBATCH --cpus-per-task=4
#SBATCH --array=1-576
#SBATCH --output=slurm_files/alignment_%a.txt

ml restore bioinformatics

bn=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/groups/dsfisher/adityam/sherlock_lab/cerevisiae/genomes/wgs_filenames.txt)

ref_genome=$"/home/groups/dsfisher/adityam/sherlock_lab/cerevisiae/S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa"
ref_genome_fasta=$"/home/groups/dsfisher/adityam/sherlock_lab/cerevisiae/S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fasta"
dir=$"/scratch/users/adityam2/wgs_data/"
out_dir=$"/scratch/users/adityam2/wgs_data/bam/"
out_dir_vcf=$"/scratch/users/adityam2/wgs_data/vcf/"

### generates bam files; puts output files in the same directory regardless of which sequencing run they were from
bwa mem -t 4 "$ref_genome" "$dir$bn"_R1_001.fastq.gz "$dir$bn"_R2_001.fastq.gz | samtools view -S -b - | samtools sort > "$out_dir${bn##*/}".sorted.bam

### runs commands to add read group information and mark duplicates
java -jar /home/groups/dsfisher/adityam/sherlock_lab/software/picard.jar AddOrReplaceReadGroups I="$out_dir${bn##*/}".sorted.bam O="$out_dir${bn##*/}".sorted.rg.bam RGID="${bn##*/}" RGLB=lib"${bn##*/}" RGPL=illumina RGPU=unit"${bn##*/}" RGSM="${bn##*/}"
samtools index "$out_dir${bn##*/}".sorted.rg.bam

java -jar /home/groups/dsfisher/adityam/sherlock_lab/software/picard.jar MarkDuplicates I="$out_dir${bn##*/}".sorted.rg.bam O="$out_dir${bn##*/}".sorted.rg.md.bam M="$out_dir${bn##*/}".md.metrics.txt
samtools index "$out_dir${bn##*/}".sorted.rg.md.bam

gatk --java-options "-Xmx4g" HaplotypeCaller -R "$ref_genome_fasta" -I "$out_dir${bn##*/}".sorted.rg.md.bam -ERC GVCF -O "$out_dir_vcf${bn##*/}".g.vcf.gz
