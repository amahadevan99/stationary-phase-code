"""
Generates bamsnap images for visual validation of candidate SNPs.

Usage:
    python gen_snaps.py <snp_coords_csv> <output_dir> [num_snaps]
    e.g. python gen_snaps.py snp_coords.csv snapshots/ 100

Inputs:
    - <snp_coords_csv>  — CSV produced by basic_filter.py (columns: file_name, location, ...)
    - BAM files in bam_path (sorted, read-group-added, duplicate-marked: *.sorted.rg.md.bam)
    - S288C R64-3-1 reference genome (.fasta)

Outputs:
    - {output_dir}/{n}.{chrom}.{pos}.{sample_abbrev}.png  — one PNG per SNP per sample,
      showing read alignments, coverage, and base identities at the variant site

Notes:
    Each image includes a moderate-coverage ancestral control sample for comparison.
    Requires bamsnap to be installed and available on PATH.

Dependencies:
    bamsnap, numpy, os, sys
"""
import os
import numpy as np
import sys

# generate snapshots of regions with SNPs, using bamsnap (requires installation)

# input csv file containing SNP coordinates to plot, e.g. snp_coords.csv
infile = sys.argv[1]
# output directory, e.g. snapshots
outdir = sys.argv[2]

# number of snapshots
if len(sys.argv)>3:
    num_snaps = int(sys.argv[3])
else:
    num_snaps = np.inf


ref_path = '/home/groups/dsfisher/adityam/sherlock_lab/cerevisiae/S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fasta'
bam_path = '/scratch/users/adityam2/wgs_data/bam/'

# plot a moderately-high coverage ancestor below each putative SNP, to act as a control
control_samps = ['20124FL-22-01-96_S97_L005','20124FL-22-01-94_S95_L005']

with open(infile,'r') as file:
    c = 0
    next(file)
    for line in file:
        split_line = line.strip().split(',')
        name,coord = split_line[0],split_line[1]
        bam_string = f'{bam_path}{name}.sorted.rg.md.bam '
        # bam_string += ' '.join([f'{bam_path}{s}.sorted.rg.md.bam' for s in control_samps])
        if control_samps[0] not in bam_string:
            bam_string += f'{bam_path}{control_samps[0]}.sorted.rg.md.bam'
        else:
            bam_string += f'{bam_path}{control_samps[1]}.sorted.rg.md.bam'

        # avoid SNPs close to the beginning of a chromosome
        # if int(coord.split(':')[1])<10000:
        #     print(coord)
        #     continue

        c += 1
        if c>num_snaps:
            break
        print(c)

        abbrev = name[name.find('-')+1:name.find('_')]
        outfile = str(c)+f'.{coord}.{abbrev}'.replace(':','.')

        # if 'chrVIII' not in coord:
        #     continue

        command = f'bamsnap -bam {bam_string} -pos {coord} -out {outdir}/{outfile}.png \
            -ref {ref_path} -draw coordinates bamplot coverage base'# -show_soft_clipped'
        os.system(command)
        os.system(f'rm {outdir}/{outfile}.png_bamsnap.log')
