"""
Maps WGS-sequenced clones to their evolution-experiment barcodes and assigns inferred fitness values.

Usage:
    Run directly as a script (no CLI arguments); configure date, directory paths,
    and experiment lists at the top of the file before running.

Inputs:
    - Paired-end FASTQ.gz from clone barcode sequencing (one pool per plate)
    - clone_info/forward_primers.csv, clone_info/reverse_primers.csv  — primer tag sequences
    - clone_info/{date}_clone_info.csv                                — clone metadata table
    - barcode_cluster_dir/{exp_id}_cluster.csv                        — bartender cluster output
    - fitmut_dir/{exp_id}_MutSeq_Result.csv                           — FitMut2 fitness results
    - fitmut_input_dir/{exp_id}_timepoints.csv                        — timepoint/cell-number params

Outputs:
    - {date}/barcodes.csv          — per-clone majority barcode(s) with read counts
    - {date}/labeled_barcodes.csv  — barcodes matched to cluster + inferred fitness
    - {date}/unique_barcodes.csv   — deduplicated version of labeled_barcodes.csv

Dependencies:
    numpy, pandas, Levenshtein (python-Levenshtein)
"""
import os,re,gzip
import numpy as np
import pandas as pd
import Levenshtein as Lev

############################################
def string_sim(s,t):
	return sum([a==b for a,b in zip(s,t)])

# sorts barcodes from list of barcodes by abundance
# also return the number of reads of each barcode, and the total number of reads associated with the primer pair
# the fitnesses used are those from the fitness inference with a higher threshold in the odds ratio
# for calling a lineage adaptive
# biggest one first
def get_sorted_barcode(bc_list,num):
	values,counts = np.unique(bc_list,return_counts=True)

	# merge barcodes similar to most abundant
	max_idx = np.argmax(counts)
	big_clust=[]
	for j,v in enumerate(values):
		if string_sim(v,values[max_idx])>len(v)-3:
			big_clust.append(j)
	v_merged = [values[j] for j in range(len(values)) if j not in big_clust]
	c_merged = [counts[j] for j in range(len(values)) if j not in big_clust]
	v_merged.append(values[max_idx])
	c_merged.append(sum(counts[big_clust]))
	v_merged = np.array(v_merged)
	c_merged = np.array(c_merged)

	idx = np.argsort(c_merged)[-num:]
	return v_merged[idx[::-1]],c_merged[idx[::-1]],np.sum(c_merged)


# the goal is to associate each primer pair (and thereby colony) with a barcode and fitness
# date = '04-24-23'
date = '08-18-23'
output_dir = date+'/'

raw_data_dir = '/scratch/users/adityam2/%s_cerevisiae/clone_barcodes/Fastq/'%date
infile_R1 = 'JT2_S1_L001_R1_001.fastq.gz'
infile_R2 = 'JT2_S1_L001_R2_001.fastq.gz'

barcode_cluster_dir = './../barcode_counts/combined_bc1_FR/'
fitmut_dir = './../fitness_inference/fitmut_output/direct_search/'
fitmut_input_dir = './../fitness_inference/fitmut_input/params0/'

############################################
forward_primer_dict = {} # maps multitag sequence to name (eg ACTGCTCT to F1)
with open('./clone_info/forward_primers.csv','r') as file:
	for line in file:
		parts = line.strip().split(',')
		forward_primer_dict[parts[1][33:41]] = parts[0]

reverse_primer_dict = {} # maps multitag sequence to name (eg ACTGCTCT to R1)
with open('./clone_info/reverse_primers.csv','r') as file:
	for line in file:
		parts = line.strip().split(',')
		reverse_primer_dict[parts[1][33:41]] = parts[0]

primer_pair_list = [] # list of primer pairs in order
exp_name_dict = {} # dict of experiment names indexed by primer pairs
with open('./clone_info/%s_clone_info.csv'%date,'r') as file:
	kk = 0
	for line in file:
		if kk==0:
			kk+=1
			continue
		parts = line.strip().split(',')
		primer_pair_list.append(parts[11])
		exp_name_dict[parts[11]] = parts[1].split(' ')[0]
		kk+=1

# writes top q majority barcode to file, along with experimental condition and primer pair
# if fewer than q barcodes for a particular primer pair, writes NaN
q = 2

############################################
# multitag_file = open('parsed.txt','w')

count = np.inf
ii=0

# looks for exact match with primer sequence
# prints the median of all barcodes which have a certain primer pair to a file

barcode_dict = {} # maps primer pair to list of barcode sequences

with gzip.open(raw_data_dir+infile_R1,'r') as gzfile1:
	with gzip.open(raw_data_dir+infile_R2,'r') as gzfile2:
		while True:
			read1 = [gzfile1.readline() for kk in range(4)]
			read2 = [gzfile2.readline() for kk in range(4)]
			fwd_line = read1[1]
			rev_line = read2[1]
			if fwd_line == b'' or ii>count:
				break

			fwd_tag = str(fwd_line[:8])[2:-1]
			rev_tag = str(rev_line[:8])[2:-1]

			# bc1 = re.findall(rb'TACC[ATCGN]{4,6}AA[ATCGN]{4,6}AA[ATCGN]{4,6}TT[ATCGN]{4,6}ATAA',fwd_line[50:100])
			# if len(bc1)>0:
			# 	bc1 = str(bc1[0][4:-4])[2:-1]

			str_fwd = str(fwd_line[:100])[2:-1]
			a = str_fwd.find('TCGGTACC')
			if a==-1:continue

			bc1 = str_fwd[a+8:a+8+26]
			# bc2 = rev_line[45-10:79+2]

			fwd_tag_id = rev_tag_id = '***'
			if fwd_tag in forward_primer_dict:
				fwd_tag_id = forward_primer_dict[fwd_tag]
			if rev_tag in reverse_primer_dict:
				rev_tag_id = reverse_primer_dict[rev_tag]
			# multitag_file.write(','.join([fwd_tag_id,rev_tag_id,bc1])+'\n')
			pp_key = rev_tag_id+fwd_tag_id
			if pp_key not in barcode_dict:
				barcode_dict[pp_key] = [bc1]
			else:
				barcode_dict[pp_key].append(bc1)
			ii += 1
# multitag_file.close()


#######################################
# writes top q majority barcode to file, along with experimental condition and primer pair
# if no barcode was seen with a particular primer pair, writes '***'

barcode_outfile = open(output_dir+'barcodes.csv','w')
barcode_outfile.write('experiment,primer pair,')
for kk in range(1,q+1):
	barcode_outfile.write('maj%s,'%kk)
for kk in range(1,q+1):
	barcode_outfile.write('nr%s,'%kk)
barcode_outfile.write('primer pair reads\n')

for pp in primer_pair_list:
	if pp in barcode_dict:
		bcs, nreads, pp_num = get_sorted_barcode(barcode_dict[pp],q)
		maj_bc = []
		for kk in range(q):
			try:
				maj_bc.append(str(bcs[kk]))
			except IndexError:
				maj_bc.append('NaN')
		for kk in range(q):
			try:
				maj_bc.append(str(nreads[kk]))
			except IndexError:
				maj_bc.append('NaN')
		maj_bc.append(str(pp_num))

		barcode_outfile.write(','.join([exp_name_dict[pp],pp,*maj_bc])+'\n')
	else:
		barcode_outfile.write(','.join([exp_name_dict[pp],pp,'***\n']))

barcode_outfile.close()
#######################################
# goes through barcode file and matches with clustered barcodes
# checks for barcode match before checking for sufficient majority

exp_name_list = ['2_day_r2','4_day_r1','4_day_r2','6_day_r1','6_day_r3','8_day_r2',
                 '8_day_r3b','10_day_r1','10_day_r2']

cluster_dict = {}
fitness_dict = {}
mut_num_dict = {}
fitmut_params = {}
for exp_id in exp_name_list:
	cluster_dict[exp_id] = pd.read_csv(barcode_cluster_dir+'%s_cluster.csv'%exp_id)
	fitness_dict[exp_id] = pd.read_csv(fitmut_dir + '%s_MutSeq_Result.csv'%exp_id)
	mut_num_dict[exp_id] = pd.read_csv(fitmut_dir + '%s_Cell_Number_Mutant_Estimated.csv'%exp_id,header=None)
	fitmut_params[exp_id] = pd.read_csv(fitmut_input_dir+'%s_timepoints.csv'%exp_id,header=None)

barcode_seqs = pd.read_csv(output_dir+'barcodes.csv')
clone_info = pd.read_csv('clone_info/%s_clone_info.csv'%date,keep_default_na=False)

outfile_labeled_bc = open(output_dir+'labeled_barcodes.csv','w')
outfile_labeled_bc.write(','.join(list(clone_info)))
outfile_labeled_bc.write(',inferred fitness,inferred mutant fraction of lineage\n')
eid=''
maj_thresh = 0.5
for index,row in barcode_seqs.iterrows():
	clone_row = [str(x) for x in np.array(clone_info.iloc[index,:])]
	if row[2]=='***':
		outfile_labeled_bc.write(','.join(clone_row[:10]))
		outfile_labeled_bc.write(',no primer pair match found\n')
		continue

	if row[0]!=eid: # takes advantage of the fact that csv is roughly organized by experiment index
		eid = row[0]
		fitmut = fitness_dict[eid]
		cluster = cluster_dict[eid]
		mut_num = mut_num_dict[eid]
		fm_params = fitmut_params[eid]
		clustered_list = np.array(cluster['Center'])

	# lev_dists = np.array([Lev.distance(xx,row[2]) for xx in clustered_list])
	# match_idxs = np.where(lev_dists<5)[0]
	match_idxs = np.nonzero(row[2]==clustered_list)[0]

	# if the majority barcode is not found in full sequencing run data
	if len(match_idxs)==0:
		outfile_labeled_bc.write(','.join(clone_row[:10]))
		outfile_labeled_bc.write(',no barcode match found\n')
		continue

	# only writes barcodes which are a substantial majority and whose runner up is sufficiently small
	if row[2+q]/row[2+2*q]<maj_thresh or row[3+q]/row[2+q]>.3:
		outfile_labeled_bc.write(','.join(clone_row[:10]))
		w1 = np.round(row[2+q]/row[2+2*q],2)
		w2 = np.round(row[3+q]/row[2+q],2)
		outfile_labeled_bc.write(',majority is %s, runner up ratio is %s, thresh is %s \n'%(w1,w2,maj_thresh))
		continue

	# assumes 8 generations per time point (but can be changed)
	mi = match_idxs[0]
	timepoint = int(clone_row[1].split(' ')[-1][1:]) # the timepoint from which the clone is sampled
	try:
		time_idx = np.where(fm_params.iloc[:,0]//8==timepoint)[0][0] # the index corresponding to this timepoint
	except IndexError: # look for next timepoint if can't find required
		try:
			time_idx = np.where(fm_params.iloc[:,0]//8==timepoint-1)[0][0]
		except IndexError:
			time_idx = np.where(fm_params.iloc[:,0]//8==timepoint+1)[0][0]

	eff_cell_num = fm_params.iloc[time_idx,1] # lineage size (effective cell number)
	num_reads = np.sum(cluster.iloc[:,time_idx+3]) # number of reads in lineage
	if cluster.iloc[mi,time_idx+3]==0:
		mut_frac = 0 # if there are no barcodes seen from the barseq at the expected timepoint
	else:
		mut_frac = np.round(mut_num.iloc[mi,time_idx]/eff_cell_num/(cluster.iloc[mi,time_idx+3]/num_reads),2)

	outfile_labeled_bc.write(','.join([*clone_row[:10],row[2],row[1],str(fitmut['Fitness'][mi]),str(mut_frac)])+'\n')

outfile_labeled_bc.close()

#####################################
# makes a separate line for each unique barcode

labeled_bc = pd.read_csv(output_dir+'labeled_barcodes.csv',keep_default_na=False)
counts,idxs = np.unique(np.array(labeled_bc['BC sequence']),return_index=True)


outfile_unique_bc = open(output_dir+'unique_barcodes.csv','w')
outfile_unique_bc.write(','.join(list(labeled_bc))+'\n')

for kk in sorted(idxs):
	if ' ' in labeled_bc.iloc[kk,10]:
		continue
	outfile_unique_bc.write(','.join([str(y) for y in labeled_bc.iloc[kk,:]])+'\n')
outfile_unique_bc.close()


