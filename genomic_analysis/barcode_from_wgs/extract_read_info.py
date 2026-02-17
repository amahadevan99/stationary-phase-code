"""
Extracts the barcode sequence from WGS BAM files via majority vote across aligned reads (QC).

Usage:
    Run directly as a script (no CLI arguments); configure bam_path at the top of the file.
    Run after align_bc.sh has produced BAM files aligned to the barcode locus.

Inputs:
    - bc_bam/*.sorted.bam  — BAM files aligned to the barcode insertion locus (from align_bc.sh)

Outputs:
    - retrieved_barcodes.csv  — two-column CSV: sample name, majority-vote barcode sequence
      '**************************' indicates no reads aligned to the barcode region
      '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' indicates reads aligned outside the barcode coordinates

Notes:
    Extracts 26 bp from positions 37-63 of the YBR209W_locus reference (the barcode coordinates).
    Soft-clipped bases are included in the majority vote to handle reads spanning the locus edge.

Dependencies:
    numpy, pysam, os, collections
"""
import numpy as np
import pysam
import os
from collections import Counter


def majority_vote(array):
    return np.array([
        Counter([base for base in column if base != '-']).most_common(1)[0][0]
        if any(base != '-' for base in column) else '^'  # Handle cases where only '-' is present
        for column in array.T  # Transpose to iterate over columns
    ])

bam_path = '/scratch/users/adityam2/wgs_data/bc_bam/'
bam_list = sorted(os.listdir(bam_path))


# ************** indicates that there were no reads that aligned to the reference in this region
# ^^^^^^^^^^^^^^ indicates that there were reads that aligned in this region, but outside the area where the barcode is at

with open('retrieved_barcodes.csv','w') as file:
    for bam in bam_list:
        if bam[-4:] == '.bai':
            continue

        # Open the BAM file
        bamfile = pysam.AlignmentFile(bam_path+bam, "rb")

        # Define the reference region
        chrom = "YBR209W_locus"  # Change to your reference sequence name
        start = 0    # Start position
        end = 70    # End position

        # Get reads in the specified region
        reads = list(bamfile.fetch(chrom, start, end))

        # Determine the length of the region
        region_length = end - start + 1

        # Create an empty list to store reads as rows
        aligned_reads = []

        # Process each read
        for read in reads:
            read_start = read.reference_start - start  # Offset to align with the region
            read_seq = read.query_sequence
            qual_seq = read.query_qualities
            cigar = read.cigar  # Get CIGAR operations

            if not cigar:
                continue
            # Initialize adjusted sequence
            adjusted_seq = list(read_seq)

            # replace low quality bases with 'Q'
            # adjusted_seq_tmp = list(read_seq)
            # adjusted_seq = []
            # for i in range(len(qual_seq)):
            #     if qual_seq[i]>30:
            #         adjusted_seq.append(adjusted_seq_tmp[i])
            #     else:
            #         adjusted_seq.append('Q')

            # Check for soft clipping (S) in the CIGAR string
            if cigar[0][0] == 4:  # Soft clipping at the start
                soft_clip_len = cigar[0][1]
                soft_clipped_bases = read_seq[:soft_clip_len]
                read_start -= soft_clip_len  # Shift read start to include soft-clipped bases
                adjusted_seq = list(soft_clipped_bases) + adjusted_seq[soft_clip_len:]  # Include clipped bases

            if cigar[-1][0] == 4:  # Soft clipping at the end
                soft_clip_len = cigar[-1][1]
                soft_clipped_bases = read_seq[-soft_clip_len:]
                adjusted_seq += list(soft_clipped_bases)  # Append clipped bases

            # Ensure read_start does not go negative
            if read_start < 0:
                adjusted_seq = adjusted_seq[abs(read_start):]  # Trim excess soft-clipped bases
                read_start = 0

            # Create an empty row filled with '-'
            row = ["-"] * region_length

            # Insert adjusted read sequence at the correct position
            sec_len = min(read_start + len(adjusted_seq),region_length) - read_start
            row[read_start:read_start+sec_len] = adjusted_seq[:sec_len]

            # Append to the list
            aligned_reads.append(row)

        # Convert to a NumPy array
        aligned_array = np.array(aligned_reads)

        if len(aligned_array)==0:
            barcode = ''.join(['*']*26)
        else:
            barcode = ''.join(majority_vote(aligned_array[:,37:63]))
        write_name = bam.split('.')[0]
        file.write(','.join([write_name,barcode])+'\n')
        print(write_name)


# save the array
# np.savez('alignment',data=aligned_array)
