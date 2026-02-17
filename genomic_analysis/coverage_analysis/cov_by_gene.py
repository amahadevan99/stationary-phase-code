"""
Computes mean sequencing coverage per annotated gene across all WGS samples.

Usage:
    Run directly as a script (no CLI arguments); configure directory paths at the top of the file.

Inputs:
    - wgs_filenames.txt                   — one sample filename prefix per line
    - coverage/{sample}.txt              — per-position coverage from gen_coverage.py
    - S288C R64-3-1 GFF annotation file  — for gene coordinates

Outputs:
    - cov_by_gene_data.npz  — numpy archive containing cov_dict:
        {sample_name: {gene_id: {len, start, sum_cov, chrom}}}
      where sum_cov / len gives mean coverage per gene per sample

Dependencies:
    numpy, pandas, biopython (Bio.SeqIO), cyvcf2, scipy
"""
import cyvcf2
import numpy as np
import sys
import gzip
import pandas as pd
import csv
from scipy.stats import binned_statistic
from Bio import SeqIO

# 1. Parse the GFF file and extract gene information
def parse_gff(gff_file):
    genes = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue  # Skip comments
            columns = line.strip().split("\t")
            chrom = columns[0]  # Chromosome name
            feature_type = columns[2]  # Feature type (e.g., gene)
            start = int(columns[3])  # Start position (1-based)
            end = int(columns[4])  # End position (1-based)
            attributes = columns[8]  # Attributes (contains gene name, etc.)

            # Extract gene ID or name from the attributes column
            gene_id = None
            for attr in attributes.split(";"):
                if attr.startswith("gene="):  # Example for GFF3
                    gene_id = attr.split("=")[1]

            if feature_type == "gene" and gene_id:
                genes.append({
                    'chrom': chrom,
                    'start': start - 1,  # Convert to 0-based indexing
                    'end': end,
                    'gene_id': gene_id
                })
    return genes


chrom_sizes = {
    "chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933,
    "chrV": 576874, "chrVI": 270161, "chrVII": 1090940, "chrVIII": 562643,
    "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177,
    "chrXIII": 924431, "chrXIV": 784333, "chrXV": 1091291, "chrXVI": 948066,"chrmt":85779
}
cumulative_offsets = {}
current_offset = 0

for chrom, size in chrom_sizes.items():
    cumulative_offsets[chrom] = current_offset
    current_offset += size

def get_genomic_position(chrom, pos):
    """Convert (chromosome, position) to a single genomic coordinate."""
    return cumulative_offsets[chrom] + pos

###########################################################################################
###########################################################################################
###########################################################################################

ref_dir = '/home/groups/dsfisher/adityam/sherlock_lab/cerevisiae/S288C_reference_genome_R64-3-1_20210421/'
gff_file = ref_dir+"saccharomyces_cerevisiae_R64-3-1_20210421.gff"

# Parse the GFF file to get genes and their coordinates
genes = parse_gff(gff_file)


# takes coverage files produced by samtools depth and bins the coverages by gene; takes average across gene

ref_genome=ref_dir+'S288C_reference_sequence_R64-3-1_20210421.fsa'
cov_dir = '/scratch/users/adityam2/wgs_data/coverage/'

with open('/home/groups/dsfisher/adityam/sherlock_lab/cerevisiae/genomes/wgs_filenames.txt', 'r') as file:
    lines = file.readlines()

cov_dict = {}
for fname in lines:
    name = fname.split('/')[1].strip()
    cov_dict[name] = {}
    print(name,flush=True)

    coverage_file = f'{cov_dir+name}.txt'
    coverage_df = pd.read_csv(coverage_file, sep="\t", header=None)
    coverage_df.columns = ["chrom", "pos", "cov"]


    for gene in genes:
        genome_start = get_genomic_position(gene["chrom"],gene["start"]) - 1
        genome_end = get_genomic_position(gene["chrom"],gene["end"])
        subset = coverage_df[genome_start:genome_end]

        # make sure that indexing is working: accessing the coverage file by cumsum index of the gene in question is valid
        assert(subset.iloc[0]['pos']==gene['start'])

        sum_coverage = subset["cov"].sum()
        gene_len = gene["end"]-gene["start"]
        cov_dict[name][gene['gene_id']] = {'len':gene_len,'start':gene["start"],'sum_cov':sum_coverage,'chrom':gene['chrom']}

np.savez('cov_by_gene_data',cov_dict=cov_dict)
