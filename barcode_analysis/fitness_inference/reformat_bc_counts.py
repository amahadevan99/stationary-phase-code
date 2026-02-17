"""
Reformats bartender cluster output into count matrices suitable for FitMut2 (or FitSeq2) input.

Usage:
    Run directly as a script (no CLI arguments); configure exp_name_list and file paths
    at the top of the file before running.

Inputs:
    - combined_bc1_FR/{exp_name}_cluster.csv  — bartender combined cluster output (one per experiment)

Outputs:
    - fitmut_input/filtered/{exp_name}_counts.csv  — barcode read-count matrix (rows = lineages,
      columns = timepoints), with low-read lineages filtered out by inverse participation ratio (IPR)

Notes:
    Lineages are filtered if their IPR across timepoints falls below ipr_thres (default 1.3).
    IPR < 1 indicates reads concentrated in a single timepoint; higher values indicate spread.

Dependencies:
    numpy, os
"""
import numpy as np
import os

# filters trajectories out if their inverse participation ratio is too small
filter_ipr = True
ipr_thres = 1.3
def ipr(x):
	s = sum(x)
	return 1/sum((x/s)**2)

exp_name_list = ['2_day_r2','4_day_r1','4_day_r2','6_day_r1','6_day_r3',
                 '8_day_r2','8_day_r3b','10_day_r1','10_day_r2']
# exp_name_list = ['8_day_r2','8_day_r3b']
for exp_name in exp_name_list:
	infile_csv = './../barcode_counts/combined_bc1_FR/%s_cluster.csv'%(exp_name)
	outfile = './fitmut_input/filtered/%s_counts.csv'%exp_name

	nlines = 0
	nfiltered = 0
	with open(outfile,'w') as outfile:
		with open(infile_csv,'r') as infile:
			for line in infile:
				if line[:7] == 'Cluster': continue
				if filter_ipr and ipr(np.array([float(x) for x in line.split(',')[3:]]))<ipr_thres:
					nfiltered += 1
					continue
				outfile.write(','.join(line.split(',')[3:]))
				nlines += 1
	print('filtered %s fraction of lineages'%(nfiltered/(nfiltered+nlines)))
