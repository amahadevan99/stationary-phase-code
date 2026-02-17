"""
Demultiplexes paired-end barcode sequencing reads, extracts barcodes, and clusters with bartender.

Usage:
    Run directly as a script (no CLI arguments); configure date_list, primer pair lists,
    and directory paths at the top of the file before submission to the Sherlock SLURM cluster.

Inputs:
    - Raw paired-end FASTQ.gz files from barcode sequencing
    - barcode_info/{date}/multitags_list.csv  — maps multitag sequences to primer names
    - barcode_info/{date}/filename_key.csv    — maps filenames to N,S primer combinations
    - barcode_info/{date}/bc_indices.csv      — maps timepoint labels to primer combinations

Outputs:
    - Demultiplexed FASTQ files split by forward/reverse primer pair
    - Extracted barcode text files (one barcode + UMI per line)
    - bartender cluster and quality CSV files per sample/primer-pair combination
    - Combined bartender output across timepoints per experiment (combined_bc1_find_loc/)

Dependencies:
    bartender (bartender_extractor_com, bartender_single_com, bartender_combiner_com), numpy
"""
import os,re,gzip
import numpy as np

date_list = ['10-06-22','02-28-23']

use_umi = True
if use_umi:
    umi_string = ''
else:
    umi_string = '_no_UMI'

# pairs of forward and reverse primer indices
p1 = ['F201','F202','F203','F204','F201','F203','F206','F207']
p2 = ['R301','R302','R303','R304','R304','R304','R304','R304']
Np = len(p1)

for date in date_list:
    #####################################
    # prepare directories
    with open('%s_barcode_files.txt'%date,'r') as file_names:
        file_list = file_names.read().split()

    # find loc tries to find the barcode based on its preceding region, while fixed loc just
    # excises a predetermined portion of the read
    raw_data_dir = '/scratch/users/adityam2/%s_cerevisiae/raw_data/'%date
    demux_dir = '/scratch/users/adityam2/%s_cerevisiae/glyc_eth/barcodes/demultiplexed_reads_FR/'%date
    extract_dir = '/scratch/users/adityam2/%s_cerevisiae/glyc_eth/barcodes/bartender/extracted_bc1_find_loc%s/'%(date,umi_string)
    cluster_dir = '/scratch/users/adityam2/%s_cerevisiae/glyc_eth/barcodes/bartender/clustered_bc1_find_loc%s/'%(date,umi_string)

    # raw_data_dir = './subset_trial/data_subset/'
    # demux_dir = './subset_trial/demultiplexed_reads/'
    # extract_dir = './subset_trial/manually_extracted/'

    for dirname in [demux_dir,extract_dir,cluster_dir]:
        if not os.path.exists(dirname):
            os.makedirs(dirname)

    #############################################
    # load information about barcode tags
    multitag_dict = {} # maps multitag sequence to name (eg ACTGCT to F201)
    with open('./barcode_info/%s/multitags_list.csv'%date,'r') as file:
        for line in file:
            parts = line.strip().split(',')
            multitag_dict[parts[1]] = parts[0]

    filename_dict = {} # maps filename to N and S names eg N701_S502
    with open('./barcode_info/%s/filename_key.csv'%date,'r') as file:
        cc = 0
        for line in file:
            if cc>0:
                parts = line.strip().split(',')
                filename_dict[parts[0]] = parts[1]+'_'+parts[3]
            cc+=1

    ###########################################
    # number of samples for which to run processing
    start_idx = 0
    end_idx = len(file_list)//2

    # prefixes of different samples (half the number of fastq files, since R1 and R2 not differentiated)
    prefix_list = []
    for idx in range(0,2*end_idx,2):
        prefix_list.append(file_list[idx][:25])

    #############################
    # the N,S primers correspond to the second step PCR primers and have already been demultiplexed into different files by Illumina.
    # Here we demultiplex reads according to first step PCR primers F,R. Write reads to different files according to the primer identity.

    # for idx in range(2*start_idx,2*end_idx,2):
    #     infile_gz1 = raw_data_dir+file_list[idx]
    #     infile_gz2 = raw_data_dir+file_list[idx+1]

    #     NS_name = filename_dict[prefix_list[idx//2][:16]]

    #     # contains information about the primers recognized in each read in the original barcode data
    #     # if these primers match the first step PCR multitags, then the read is written to the files in
    #     # outfile_R1_fastq_dict and outfile_R2_fastq_dict
    #     multitag_file = open(demux_dir+NS_name+'.txt','w')

    #     outfile_R1_fastq_dict = {}
    #     outfile_R2_fastq_dict = {}

    #     for x in range(Np):
    #         outfile_R1_fastq_dict[x] = open(demux_dir+NS_name+'_%s_%s_R1.fastq'%(p1[x],p2[x]),'wb')
    #         outfile_R2_fastq_dict[x] = open(demux_dir+NS_name+'_%s_%s_R2.fastq'%(p1[x],p2[x]),'wb')

    #     with gzip.open(infile_gz1,'r') as gzfile1:
    #         with gzip.open(infile_gz2,'r') as gzfile2:
    #             while True:
    #                 read1 = [gzfile1.readline() for kk in range(4)]
    #                 read2 = [gzfile2.readline() for kk in range(4)]
    #                 fwd_line = read1[1]
    #                 rev_line = read2[1]
    #                 if fwd_line == b'':
    #                     break

    #                 fwd_tag = str(fwd_line[8:8+6])[2:-1]
    #                 rev_tag = str(rev_line[8:8+9])[2:-1]

    #                 fwd_tag_id = rev_tag_id = '****'
    #                 if fwd_tag in multitag_dict:
    #                     fwd_tag_id = multitag_dict[fwd_tag]
    #                 if rev_tag in multitag_dict:
    #                     rev_tag_id = multitag_dict[rev_tag]
    #                 multitag_file.write(fwd_tag+' '+fwd_tag_id+' ')
    #                 multitag_file.write(rev_tag+' '+rev_tag_id+'\n')

    #                 for x in range(Np):
    #                     if fwd_tag_id == p1[x] and rev_tag_id == p2[x]:
    #                         for l in read1: outfile_R1_fastq_dict[x].write(l)
    #                         for l in read2: outfile_R2_fastq_dict[x].write(l)
    #                         break
    #     multitag_file.close()
    #     for x in outfile_R1_fastq_dict:
    #         outfile_R1_fastq_dict[x].close()
    #     for x in outfile_R2_fastq_dict:
    #         outfile_R1_fastq_dict[x].close()

    #############################
    # extract barcodes and UMIs manually (can also do with bartender)
    for idx in range(start_idx,end_idx):
        for x in range(Np):
            primer_names = filename_dict[prefix_list[idx][:16]]+'_%s_%s'%(p1[x],p2[x])
            print('extracting barcodes from %s'%primer_names,flush=True)
            infile_fq1 = demux_dir+primer_names+'_R1.fastq'
            infile_fq2 = demux_dir+primer_names+'_R2.fastq'
            outfile_txt = extract_dir+primer_names
            clusterfile = cluster_dir+primer_names

            # extract forward UMI and 1st barcode with bartender
            # os.system('bartender_extractor_com -f {} -o {} -q ? -p TACC[4-7]AA[4-7]AA[4-7]TT[4-7]ATAA -m \
            #     2 --direction=forward -u 0,8'.format(infile_fq1,outfile_txt))

            # extract forward and reverse UMI and both barcodes
            # only write first barcode to file
            with open(outfile_txt+'_barcode.txt','wb') as extracted_file:
                fq1 = open(infile_fq1,'rb')
                fq2 = open(infile_fq2,'rb')
                counter = 0
                while True:
                    read1 = [fq1.readline() for kk in range(4)]
                    read2 = [fq2.readline() for kk in range(4)]
                    fwd_line = read1[1]
                    rev_line = read2[1]
                    if fwd_line == b'':
                        break
                    fwd_umi = fwd_line[:8]
                    rev_umi = rev_line[:8]

                    # extract barcode from fixed location (fixed loc)
                    # bc1 = fwd_line[59+4:93-4]

                    # extract barcode based on position of flanking sequence (find loc)
                    a = fwd_line.find(b'TATCGGTACC')
                    # if a!=55:print(a,flush=True)
                    if a==-1 or a>110:continue # avoid hitting end of read
                    bc1 = fwd_line[a+10:a+10+26]

                    if use_umi:
                        extracted_file.write(b','.join([bc1,b''.join([fwd_umi,rev_umi])])) # with UMIs
                    else:
                        extracted_file.write(b','.join([bc1,str(counter).encode()])) # no UMIs
                    extracted_file.write(b'\n')
                    counter += 1
                fq1.close()
                fq2.close()

    ##########################
    # cluster barcodes with bartender
    for idx in range(start_idx,end_idx):
        for x in range(Np):
            primer_names = filename_dict[prefix_list[idx][:16]]+'_%s_%s'%(p1[x],p2[x])
            outfile_txt = extract_dir+primer_names
            clusterfile = cluster_dir+primer_names
            print('clustering barcodes from %s'%primer_names,flush=True)
            command = 'bartender_single_com -f {} -o {} -d 3'.format(outfile_txt+'_barcode.txt',clusterfile)
            print(command)
            os.system(command)

########################################################
########################################################
# combine timepoints from the same experiment
combined_dir = './combined_bc1_find_loc%s/'%(umi_string)
if not os.path.exists(combined_dir):
    os.makedirs(combined_dir)


timepoint_dict = {} # maps description to N,S,F,R combination eg 4_day_r2T1 to N701_S502_F203_R303
directory_dict = {} # maps description to folder containing the clustered barcodes (which depends on sequencing date)
for dd in date_list:
    with open('./barcode_info/%s/bc_indices.csv'%dd,'r') as file:
        for line in file:
            comp = line.strip().split(',')
            key = comp[0].replace(' ','_')
            value = '_'.join(comp[1:])
            timepoint_dict[key] = value
            directory_dict[key] = ('/scratch/users/adityam2/%s_cerevisiae/glyc_eth'
                '/barcodes/bartender/clustered_bc1_find_loc%s/'%(dd,umi_string))

################################
exper_ids = ['2_day_r2','4_day_r1','4_day_r2','6_day_r1','6_day_r3','8_day_r2','8_day_r3b','10_day_r1','10_day_r2']
# exper_ids = ['8_day_r2','8_day_r3b']

if use_umi:
    suffixes = ['_pcr_cluster.csv','_pcr_quality.csv'] # with UMIs
else:
    suffixes = ['_cluster.csv','_quality.csv'] # no UMIs

def get_filekey(eid,i):
    if i==0:
        return eid.split('r')[0]+'T0'
    if eid=='8_day_r3b' and i<=5:
        return eid[:-1]+'T%s'%i
    else:
        return eid+'T%s'%i

for kk in range(len(exper_ids)):
    eid = exper_ids[kk]
    tp = 20 # only takes files whose keys are in timepoint_dict
    comb_file_list = []
    for i in range(2*tp):
        if get_filekey(eid,i//2) in timepoint_dict:
            fkey = get_filekey(eid,i//2)
            comb_file_list.append(directory_dict[fkey]+timepoint_dict[fkey]+suffixes[i%2])
    os.system('bartender_combiner_com -f {} -o {} -c 5'.format(','.join(comb_file_list),combined_dir+eid))
    # choosing -c 5 here means that there can be differences even if the same file is combined with other files,
    # because of the content of the other files.


