"""
Demultiplexes and extracts barcodes from fitness remeasurement FASTQ files.

Usage:
    python extract_bcs.py <date>
    e.g. python extract_bcs.py 07-17-24

Inputs:
    - {date}_barcode_files.txt                     — list of raw FASTQ.gz filenames
    - Raw paired-end FASTQ.gz files in raw_data_dir
    - barcode_info/primer_indices.csv              — FOS/ROS primer tag sequences

Outputs:
    - demux_dir/{FOS}_{ROS}_R1.txt / _R2.txt      — demultiplexed read sequences by primer pair
    - extract_dir/{FOS}_{ROS}_barcode.txt          — extracted barcode sequences (one per line)

Dependencies:
    numpy, os, sys, gzip
"""
import os,sys,gzip
import numpy as np

# MAKE SURE DEMULTIPLEX DIRECTORY is empty when running this script, otherwise it will not work
# it appends files to demultiplexed read lists, so will not overwrite existing files

date = sys.argv[1]#'07-17-24'

def hamming_dist(x,y):
    return sum(c1!=c2 for c1,c2 in zip(x,y))

# forward and reverse primer indices
p1 = ['FOS%s'%i for i in range(1,11)]
p2 = ['ROS%s'%i for i in range(1,14)]

#####################################
# prepare directories
with open(f'{date}_barcode_files.txt','r') as file_names:
    file_list = file_names.read().split()

raw_data_dir = f'/scratch/users/adityam2/{date}_cerevisiae/raw_data/'
demux_dir = f'/scratch/users/adityam2/{date}_cerevisiae/barcodes/demultiplexed_reads/'
extract_dir = f'/scratch/users/adityam2/{date}_cerevisiae/barcodes/extracted_bc/'

# raw_data_dir = './data_subset/raw_data/'
# demux_dir = './data_subset/demux_reads/'
# extract_dir = './data_subset/extracted_barcodes/'

for dirname in [demux_dir,extract_dir]:
    if not os.path.exists(dirname):
        os.makedirs(dirname)

#############################################
# load information about primer tags
FOS_dict = {} # maps multitag sequence to name (eg ACTGCT to FOS1)
ROS_dict = {} # maps multitag sequence to name (eg ACTGCT to ROS2)
with open(f'./barcode_info/primer_indices.csv','r') as file:
    for c,line in enumerate(file):
        parts = line.strip().split(',')
        if c<10:
            FOS_dict[parts[2]] = parts[0]
        else:
            ROS_dict[parts[2]] = parts[0]

###########################################
# number of samples for which to run processing
start_idx = 0
end_idx = len(file_list)//2

# prefixes of different samples (half the number of fastq files, since R1 and R2 not differentiated)
prefix_list = []
for idx in range(0,2*end_idx,2):
    prefix_list.append(file_list[idx][:25])

#############################
# Here we demultiplex reads according to the FOS and ROS primers. Write reads to different files according to the primer identity.

for idx in range(2*start_idx,2*end_idx,2):
    infile_gz1 = raw_data_dir+file_list[idx]
    infile_gz2 = raw_data_dir+file_list[idx+1]

    # contains information about the primers recognized in each read in the original barcode data
    # if these primers match the first step PCR multitags, then the read is written to the files in
    # outfile_R1_fastq_dict and outfile_R2_fastq_dict
    multitag_file = open(demux_dir+'multitags.txt','w')

    R1_list_dict = {} # list of reads for each primer combo
    R2_list_dict = {}

    nreads = 0
    with gzip.open(infile_gz1,'r') as gzfile1:
        with gzip.open(infile_gz2,'r') as gzfile2:
            while True:
                read1 = [gzfile1.readline() for kk in range(4)]
                read2 = [gzfile2.readline() for kk in range(4)]
                fwd_line = read1[1]
                rev_line = read2[1]
                if fwd_line == b'':
                    break

                fwd_tag_id = rev_tag_id = '****'
                for FOS in FOS_dict:
                    if hamming_dist(FOS,fwd_line[:len(FOS)].decode('utf-8'))<3:
                        fwd_tag_id = FOS_dict[FOS]
                        break
                for ROS in ROS_dict:
                    if hamming_dist(ROS,rev_line[:len(ROS)].decode('utf-8'))<3:
                        rev_tag_id = ROS_dict[ROS]
                        break

                multitag_file.write(fwd_tag_id+' '+rev_tag_id+'\n')
                tag_combo = fwd_tag_id+'_'+rev_tag_id

                if tag_combo not in R1_list_dict and '*' not in tag_combo:
                    R1_list_dict[tag_combo] = []
                    R2_list_dict[tag_combo] = []
                if '*' not in tag_combo:
                    R1_list_dict[tag_combo].append(fwd_line)
                    R2_list_dict[tag_combo].append(rev_line)
                nreads += 1

                # periodically write reads to files to keep low memory profile
                if nreads%1000000==0:
                    print('%s reads processed'%nreads,flush=True)
                    for k in R1_list_dict:
                        with open(demux_dir+f'{k}_R1.txt','ab') as demux_file:
                            for l in R1_list_dict[k]:
                                demux_file.write(l)
                        with open(demux_dir+f'{k}_R2.txt','ab') as demux_file:
                            for l in R2_list_dict[k]:
                                demux_file.write(l)
                    R1_list_dict = {}
                    R2_list_dict = {}
    multitag_file.close()

    for k in R1_list_dict:
        with open(demux_dir+f'{k}_R1.txt','ab') as demux_file:
            for l in R1_list_dict[k]:
                demux_file.write(l)
        with open(demux_dir+f'{k}_R2.txt','ab') as demux_file:
            for l in R2_list_dict[k]:
                demux_file.write(l)


# #############################
# extract barcodes manually

# primer combinations that have already been parsed
tag_combo_done = set()

for idx in range(start_idx,end_idx):
    fnames = sorted(os.listdir(demux_dir))
    for fname in fnames:
        if fname=='multitags.txt':
            continue

        tag_combo = fname.split('.')[0][:-3]
        if tag_combo in tag_combo_done:
            continue
        else:
            tag_combo_done.add(tag_combo)

        print(f'extracting barcodes from {tag_combo}',flush=True)
        infile_fq1 = demux_dir+tag_combo+'_R1.txt'
        infile_fq2 = demux_dir+tag_combo+'_R2.txt'
        outfile_txt = extract_dir+tag_combo

        # extract barcode
        with open(outfile_txt+'_barcode.txt','wb') as extracted_file:
            fq1 = open(infile_fq1,'rb')
            # fq2 = open(infile_fq2,'rb')
            counter = 0
            while True:
                fwd_line = fq1.readline()
                # rev_line = fq2.readline()
                if fwd_line == b'':
                    break

                bc1 = b''
                bc2 = b''
                # extract barcode based on position of flanking sequence
                a1 = fwd_line.find(b'TCGGTACC')
                if not(a1==-1 or a1>110): # avoid hitting end of read
                    bc1 = fwd_line[a1+8:a1+8+26]

                a2 = fwd_line.find(b'GAAGTTAT')
                if not(a2==-1 or a2>110): # avoid hitting end of read
                    bc2 = fwd_line[a2+8:a2+8+26]

                if not(bc1==b'' and bc2==b''):
                    extracted_file.write(b','.join([bc1,bc2,str(counter).encode()]))
                    extracted_file.write(b'\n')


                counter += 1
            fq1.close()
            # fq2.close()


