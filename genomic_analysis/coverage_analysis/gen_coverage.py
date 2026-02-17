"""
Generates per-position sequencing coverage files for each WGS sample using samtools depth.

Usage:
    python gen_coverage.py <genome_index>
    e.g. python gen_coverage.py 0
      genome_index: integer index into wgs_filenames.txt (0-based)

Inputs:
    - ../wgs_filenames.txt                 — one sample filename prefix per line
    - bam/{sample}.sorted.rg.md.bam       — duplicate-marked BAM (from alignment.sh)

Outputs:
    - coverage/{sample}.txt  — tab-separated file with columns: chrom, position, coverage depth
      (produced by `samtools depth -a`, so all positions are included even if coverage is zero)

Dependencies:
    samtools, numpy, os, sys
"""
import os
import numpy as np
import sys

genome_idx = int(sys.argv[1])

# generate coverage files using samtools

ref_genome='/home/groups/dsfisher/adityam/sherlock_lab/cerevisiae_data/\
                S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa'
bam_dir = '/scratch/users/adityam2/wgs_data/bam/'
cov_dir = '/scratch/users/adityam2/wgs_data/coverage/'

with open('./../wgs_filenames.txt', 'r') as file:
    lines = file.readlines()


fname = lines[genome_idx]
name = fname.split('/')[1].strip()
print(name)

# filter out secondary alignements, then use samtools depth to calculate coverage
# os.system(f'samtools view -b -F 0x100 {bam_dir+name}.sorted.bam > {bam_dir+name}.primary.bam')

os.system(f'samtools depth -a {bam_dir+name}.sorted.rg.md.bam > {cov_dir+name}.txt')
