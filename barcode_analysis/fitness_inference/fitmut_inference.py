"""
Runs FitMut2 fitness inference on barcode count data.

Usage:
    python fitmut_inference.py <exp_index> <filtered>
    e.g. python fitmut_inference.py 0 1
      exp_index: integer index into exp_name_list
      filtered:  0 = use unfiltered counts, 1 = use IPR-filtered counts (from reformat_bc_counts.py)

Inputs:
    - fitmut_input/{filtered_str}/{exp_name}_counts.csv    — barcode count matrix
    - fitmut_input/params0/{exp_name}_timepoints.csv       — timepoint and cell-number parameters

Outputs:
    - fitmut_output_trial/{filtered_str}/{exp_name}_MutSeq_Result.csv
    - fitmut_output_trial/{filtered_str}/{exp_name}_Cell_Number_Mutant_Estimated.csv
    (and other FitMut2 output files)

Dependencies:
    FitMut2 (fitmut2_run.py), numpy, os, sys
"""
import numpy as np
import os
import sys

exp_name_list = ['2_day_r2','4_day_r1','4_day_r2','6_day_r1','6_day_r3',
                 '8_day_r2','8_day_r3b','10_day_r1','10_day_r2']

exp_name = exp_name_list[int(sys.argv[1])]
filtered_str = ['unfiltered','filtered'][int(sys.argv[2])]

script_file = '/home/groups/dsfisher/adityam/sherlock_lab/software/FitMut2/main_code/fitmut2_run.py'
input_file = './fitmut_input/%s/%s_counts.csv'%(filtered_str,exp_name)
timepoint_file = './fitmut_input/params0/%s_timepoints.csv'%exp_name
output_file = './fitmut_output_trial/%s/%s'%(filtered_str,exp_name)

# os.system('python3 {} -i {} -t {} -o {}'.format(script_file,input_file,timepoint_file,output_file))
os.system('python3 {} -i {} -t {} -o {} --opt_algorithm direct_search -p 0'.format(
	script_file,input_file,timepoint_file,output_file))
