"""
Counts occurrences of known barcodes in extracted barcode files from fitness remeasurement.

Usage:
    python count_bcs.py <date>
    e.g. python count_bcs.py 07-17-24

Inputs:
    - barcode_info/barcode_list.csv          — reference list of expected barcode sequences
    - extracted_bc/ directory (produced by extract_bcs.py) containing *_barcode.txt files

Outputs:
    - bc_counts_{date}.csv  — matrix of barcode counts; rows are primer pairs, columns are barcodes

Dependencies:
    numpy, csv, os, sys
"""
import os,re,gzip,sys
import numpy as np
import csv

date = sys.argv[1]#'07-17-24'
# go through each file of barcodes belonging to a particular primer pair and count the number of each of the
# expected barcodes within that file

data_dir = f'/scratch/users/adityam2/{date}_cerevisiae/barcodes/extracted_bc/'
outfile = f'bc_counts_{date}.csv'

# contains a list of the barcodes that we are looking for in each file
bc_list = []

with open('barcode_info/barcode_list.csv', mode='r', newline='') as csvfile:
    csv_reader = csv.reader(csvfile)
    for row in csv_reader:
        if row and ('BC' not in row[0]):
            bc_list.append(row[2])

# print(bc_list)


file_list = os.listdir(data_dir)

with open(outfile,'w') as file:
    file.write(','+','.join(bc_list)+'\n')

for file_name in sorted(file_list):
    print('counting %s'%file_name)
    with open(data_dir+file_name, 'r') as file:
        bc_count_dict = {b:0 for b in bc_list}
        for line in file:
            split_line = line.split(',')
            bc1 = split_line[0]
            bc2 = split_line[1]

            if bc1 in bc_count_dict:
                bc_count_dict[bc1]+=1
            elif bc1[:-1] in bc_count_dict:
                bc_count_dict[bc1[:-1]]+=1
            elif bc1[1:] in bc_count_dict:
                bc_count_dict[bc1[1:]]+=1
    # print(bc_count_dict)

    count_string = ','.join([str(bc_count_dict[b]) for b in bc_list])
    with open(outfile,'a') as file:
        primer_pair = file_name.split('_barcode')[0]+','
        file.write(primer_pair+count_string+'\n')
