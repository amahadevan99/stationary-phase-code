"""
Computes windowed sequencing coverage in a specified genomic region across all WGS samples.

Usage:
    Run directly as a script (no CLI arguments); configure chrom_id, start_pos, end_pos,
    bin_size, and directory paths at the top of the file before running.

Inputs:
    - wgs_filenames.txt          — one sample filename prefix per line
    - coverage/{sample}.txt     — per-position coverage from gen_coverage.py

Outputs:
    - regional_cov.npz  — numpy archive containing:
        cov_dict: {sample_name: {chrom: {mean_coverage, bin_edges, bin_centers, bases, total_coverage}}}
        bin_size: the window size used

Notes:
    By default analyzes chrIII:200000-250000 with 100 bp bins; edit chrom_id/start_pos/end_pos
    at the top of the script to target a different region.

Dependencies:
    numpy, scipy, matplotlib, os
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic


# takes coverage files produced by samtools depth and bins the positions and coverages for each chromosome, taking means within bins

chrom_id = 'chrIII'
start_pos = 200000
end_pos = 250000
bin_size = 100

ref_genome='/home/groups/dsfisher/adityam/sherlock_lab/cerevisiae/S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa'
cov_dir = '/scratch/users/adityam2/wgs_data/coverage/'

with open('/home/groups/dsfisher/adityam/sherlock_lab/cerevisiae/genomes/wgs_filenames.txt', 'r') as file:
    lines = file.readlines()

cov_dict = {}
for fname in lines:
    name = fname.split('/')[1].strip()
    print(name,flush=True)

    # code adapted from chatGPT
    coverage_file = f'{cov_dir+name}.txt'

    # Read coverage data from file
    chromosomes = {}
    with open(coverage_file, 'r') as f:
        for line in f:
            chrom, pos, cov = line.strip().split()
            if chrom != chrom_id or not(start_pos<int(pos)<end_pos):
                continue
            if chrom not in chromosomes:
                chromosomes[chrom] = {'positions': [], 'coverage': []}
            chromosomes[chrom]['positions'].append(int(pos))
            chromosomes[chrom]['coverage'].append(int(cov))

    # for each window of size bin_size, calculates the mean coverage
    # the last window may stretch beyond the edge of the chromosome
    for chrom, data in chromosomes.items():
        bins = np.arange(np.min(data['positions']),np.max(data['positions'])+bin_size,bin_size)
        total_cov,_,_ = binned_statistic(data['positions'],data['coverage'],statistic='sum',bins=bins)
        mean_cov_binned = total_cov/bin_size

        if name not in cov_dict:
            cov_dict[name] = {}

        cov_dict[name][chrom] = {'mean_coverage':mean_cov_binned,'bin_edges':bins,'bin_centers':(bins[1:]+bins[:-1])/2,
                'bases':len(data['positions']),'total_coverage':sum(data['coverage'])}

np.savez('regional_cov',cov_dict=cov_dict,bin_size=bin_size)
