"""
Runs FitSeq2 fitness inference on barcode count data (primary inference method).

Usage:
    python fitseq_inference.py <exp_index>
    e.g. python fitseq_inference.py 0
      exp_index: integer index into exp_name_list

Inputs:
    - fitseq_input/{exp_name}_traj.csv        — barcode read-count trajectories
    - fitseq_input/{exp_name}_timepoints.csv  — timepoints in generations

Outputs:
    - fitseq_output/{exp_name}_FitSeq2_Result.csv  — inferred fitness values per lineage
    (and other FitSeq2 output files)

Notes:
    Uses -dt 8 (8 generations per transfer cycle). Runs in single-process mode (-p 0).
    Uses a modified version of FitSeq2 adapted for this experimental design.

Dependencies:
    FitSeq2 (FitSeq2.py), numpy, os, sys
"""
import numpy as np
import os
import sys

exp_name_list = ['M05_10_day_r1','M05_10_day_r2','M05_10_day_r3','M05_2_day_r1','M05_2_day_r2',
                    'M05_2_day_r3','M05_4_day_r1','M05_4_day_r2','M05_4_day_r3','M05_6_day_r1',
                    'M05_6_day_r2','M05_6_day_r3','M05_8_day_r1','M05_8_day_r2','M05_8_day_r3',
                    'M3_1_day_r1','M3_1_day_r2','M3_1_day_r3','M3_2_day_r1','M3_2_day_r2',
                    'M3_2_day_r3','M3_3_day_r1','M3_3_day_r2','M3_3_day_r3','M3_5_day_r1',
                    'M3_5_day_r2','M3_5_day_r3']

exp_name = exp_name_list[int(sys.argv[1])]

script_file = '/home/groups/dsfisher/adityam/sherlock_lab/software/FitSeq2/src/FitSeq2.py'
input_file = f'fitseq_input/{exp_name}_traj.csv'
timepoint_file = f'fitseq_input/{exp_name}_timepoints.csv'
output_file = f'fitseq_output/{exp_name}'

os.system(f'python3 {script_file} -i {input_file} -t {timepoint_file} -o {output_file} -dt 8 -p 0')
