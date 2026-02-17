"""
Concatenates FASTQ files from multiple whole-genome sequencing runs for the same clone.

Usage:
    Run directly as a script (no CLI arguments); configure directory paths and the
    resequenced_wgs.csv metadata file at the top of the file before running.

Inputs:
    - ../filename_info/resequenced_wgs.csv  — table mapping clone indices to file keys
      from two sequencing runs (columns: clone index, file_key_run1, file_key_run2)
    - raw_data_07-17-24/, raw_data_06-12-23/, raw_data_10-23-23/  — directories of FASTQ.gz files

Outputs:
    - merged_files/{filename}_R1_001.fastq.gz  — concatenated R1 reads
    - merged_files/{filename}_R2_001.fastq.gz  — concatenated R2 reads

Notes:
    Uses `cat` to concatenate gzipped files directly (preserves gzip format).
    Intended for clones that were sequenced twice to achieve higher coverage.

Dependencies:
    pandas, numpy, os
"""
import os
import pandas as pd
import numpy as np

# concatenate the .fastq.gz files arising from the same clone. Putting all the reads together allows us to
# call SNPs from a larger collection of data


reseq = pd.read_csv('../filename_info/resequenced_wgs.csv')
reseq_dir = '/scratch/users/adityam2/wgs_data/raw_data_07-17-24/'
dir1 = '/scratch/users/adityam2/wgs_data/raw_data_06-12-23/'
dir2 = '/scratch/users/adityam2/wgs_data/raw_data_10-23-23/'
fnames1 = os.listdir(dir1)
fnames2 = os.listdir(dir2)
reseq_fnames = os.listdir(reseq_dir)

c = 0
outfile_list = []

for i in range(480):
    if i%10==0:
        print(i)
    if reseq.iloc[i,2] is np.nan:
        continue
    file1_key = reseq.iloc[i,1]
    file2_key = reseq.iloc[i,2]
    if '-19-' in file1_key:
        fnames = fnames1
        file_dir = dir1
    elif '-21-' in file1_key:
        fnames = fnames2
        file_dir = dir2

    file1_name = file2_name = 0

    for f in fnames:
        if file1_key+'_' in f and '_R1_' in f:
            file1_name = f
            break

    for f in reseq_fnames:
        if file2_key+'_' in f and '_R1_' in f:
            file2_name = f
            break

    file1 = file_dir+file1_name
    file2 = reseq_dir+file2_name
    outfile = '/scratch/users/adityam2/wgs_data/merged_files/'+file1_name
    outfile_list.append(outfile)

    command = f'cat {file1} {file2} > {outfile}'
    # command = f'zcat {file1} {file2} | gzip > {outfile}'
    # print(command)
    os.system(command)
    os.system(command.replace('_R1_','_R2_'))
print(len(np.unique(outfile_list)))
