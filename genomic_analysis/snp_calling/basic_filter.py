"""
Hard-filters a joint-genotyped VCF and outputs SNP coordinates as a CSV file.

Usage:
    Run directly as a script (no CLI arguments); configure input/output filenames
    and filter thresholds at the top of the file before running.

Inputs:
    - wgs_all_samples.vcf.gz          — joint-genotyped VCF (from merge_vcf.sh / GenotypeGVCFs)
    - ../../filename_info/evo_cond_map.csv  — maps sample filenames to evolutionary condition metadata
    - ../../coverage/mean_cov.csv     — mean coverage per sample (for reporting)
    - S288C R64-3-1 GFF annotation file

Outputs:
    - ../bamsnap/snp_coords.csv       — filtered SNP coordinates with sample/gene/condition metadata
    - ../bamsnap/dup_snp_coords.csv   — SNPs called in more than one sample (for manual inspection)
    - filtered_snps.vcf.gz            — VCF containing only passing SNPs
    - filtered_out.csv                — sites removed by multiplet or depth filters

Hard filter thresholds (GATK annotations):
    QD > 4, FS < 40, SOR < 2, MQ > 50, MQRankSum > -3, ReadPosRankSum > -5

Dependencies:
    cyvcf2, numpy, pandas, csv, gzip, biopython (Bio.SeqIO)
"""
import cyvcf2
import numpy as np
import sys
import gzip
import pandas as pd
import csv


#######################################################################################################
#######################################################################################################
#######################################################################################################
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

# 2. Find the gene at the given position
def find_gene_at_position(genes, chrom, position):
    for gene in genes:
        if gene['chrom'] == chrom and gene['start'] <= position < gene['end']:
            return gene
    return None


ref_dict = '/home/groups/dsfisher/adityam/sherlock_lab/cerevisiae/S288C_reference_genome_R64-3-1_20210421/'
# genome_fasta_file = ref_dict+"S288C_reference_sequence_R64-3-1_20210421.fasta"
gff_file = ref_dict+"saccharomyces_cerevisiae_R64-3-1_20210421.gff"

# Parse the GFF file to get genes and their coordinates
genes = parse_gff(gff_file)

barcode_dict = {} # maps filename to barcode
reader = csv.DictReader(open('../../filename_info/evo_cond_map.csv'))
for row in reader:
    barcode_dict[row['origin_track']] = row['clone_num']+','+row['barcode']

cov_dict = {} # maps filename to barcode
with open('../../coverage/mean_cov.csv') as file:
    for line in file:
        fname,cov = line.split(',')
        cov_dict[fname] = str(np.round(float(cov),3))


# reading VCF
def get_vcf_names(vcf_path):
    with gzip.open(vcf_path, "rt") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                vcf_names = [x.strip('\n') for x in line.split('\t')]
                break
    ifile.close()
    return vcf_names


# returns True if the variant fails to make it through any of the filters
def hard_filter(variant,hard_filters):

    variant_annotations = {annotation:variant.INFO.get(annotation) for annotation in hard_filters.keys()}
    variant_annotations = {key:(val if val!=None else np.nan) for key,val in variant_annotations.items()}

    filter_out = np.any([variant_annotations['QD']<hard_filters['QD'],
            (variant_annotations['FS']>hard_filters['FS']),
            (variant_annotations['SOR']>hard_filters['SOR']),
            (variant_annotations['MQ']<hard_filters['MQ']),
            (variant_annotations['MQRankSum']<hard_filters['MQRankSum']),
            (variant_annotations['ReadPosRankSum']<hard_filters['ReadPosRankSum'])])

    return filter_out

#######################################################################################################
#######################################################################################################
#######################################################################################################

# input file name, e.g. wgs_all_samples_3X_coverage.vcf.gz
# this script filters this file based on specified filters, and outputs a csv file with the
# filtered mutations, as a well as a file that IGV can read to take pictures of the desired regions.
# see documention for VCF https://gemini.readthedocs.io/en/latest/content/database_schema.html

file_name = 'wgs_all_samples.vcf.gz'#sys.argv[1]
filtered_file = 'filtered_out.csv'
outfile = 'snp_coords'#sys.argv[2] # output file, as csv, in bamsnap directory
multip_outfile = 'dup_snp_coords'#sys.argv[2] # output file, as csv, in bamsnap directory
names = get_vcf_names(file_name)
sample_names = np.array(names[9:])
evo_cond = pd.read_csv('../../filename_info/evo_cond_map.csv')


# set filter levels
QD_thresh = 4 # get rid of less than; GATK default is 2 (value has big effect on number of filtered mutations)
FS_thresh = 40 # get rid of greater than; GATK default is 60 (value has no effect on number of filtered mutations)
SOR_thresh = 2 # get rid of greater than; GATK default is 3
MQ_thresh = 50 # get rid of less than; GATK default is 40 (value has big effect on number of filtered mutations)
MQRankSum_thresh = -3 # get rid of less than; GATK default is -12.5
ReadPosRankSum_thresh = -5 # get rid of less than; GATK default is -8.0

# put them in a dictionary that will be useful later
hard_filters = {'QD':QD_thresh, # quality by depth
                'FS':FS_thresh, # fisher strand
                'SOR':SOR_thresh, # strand odds ratio
                'MQ':MQ_thresh, # mapping quality
                'MQRankSum':MQRankSum_thresh,
                'ReadPosRankSum':ReadPosRankSum_thresh
               }


# There is a hard filtering step
# there is also an ancestral filtering step, though this was maybe filtering out real mutations and can be commented out
# then there is a step that filters out things with insufficient variant depth (<2), and things which occur in more than one sample
# possibly we should keep things that occur in only e.g. 2 samples as there can be two different mutations
# in different samples at the same site

# hard filtering
vcf = cyvcf2.VCF(file_name)
variable_sites = []
hard_filtered = []
multiplet_filtered = []
insufficient_var_depth = []
ancestral_variants = []

GQ_thresh = 30 # determines which mutations are filtered out due to low quality in variant.gt_quals

num_var = 0
for variant in vcf:
    num_var += 1
    if hard_filter(variant,hard_filters):
        hard_filtered.append(variant)
        continue

    high_quality_indices = np.where(variant.gt_quals > GQ_thresh)[0]

    # if variant.end==208458:
    #     print(variant.gt_types[variant.gt_types!=0])
    #     print(sample_names[variant.gt_types!=0])
    #     print(variant.gt_quals[variant.gt_types!=0])

    ## putative ancestral mutation filtering
    # if len(high_quality_indices)>5 and len(set(variant.gt_types[high_quality_indices])) == 1:
    #     ancestral_variants.append(variant)
    #     continue

    variable_sites.append(variant)

print('initial number of variants: %s'%num_var)


#########################################################################
vcf_reader = cyvcf2.VCF('wgs_all_samples.vcf.gz')
vcf_writer = cyvcf2.Writer('filtered_snps.vcf.gz', vcf_reader)
snp_count = 0
# write the called SNPs to a file, excluding SNPs which were called in too many genomes
# additionally write mutations found in more than one sample to a separate file for special inspection
with open('../bamsnap/%s.csv'%outfile,'w') as file:
    with open('../bamsnap/%s.csv'%multip_outfile,'w') as multip_file:
        file.write('file_name,location,gene,evol_condition,clone_num,barcode,mean_cov,pos\n')
        multip_file.write('file_name,location,gene,evol_condition,clone_num,barcode,mean_cov,pos\n')
        for variant in variable_sites:
            # the stipulation that the variant depth is bigger than 2 filters out a lot
            alt_locs = np.where((variant.gt_types!=0) * (variant.gt_alt_depths>2))[0]
            samp_list = sample_names[alt_locs]

            # filter out mutations which don't appear with sufficient depth
            if len(samp_list)==0:
                insufficient_var_depth.append(variant)
                continue

            # filter out mutations which are called in more than 3 samples
            if len(samp_list)>3:
                multiplet_filtered.append((variant,samp_list))
                continue

            # filter out mutations where the same allele is seen in all samples
            if len(samp_list)>1:
                mut_alleles = np.array([variant.gt_bases[s] for s in alt_locs])
                if len(np.unique(mut_alleles)) == 1:
                    multiplet_filtered.append((variant,samp_list))
                    continue

            for samp_string in samp_list:
                if '-22-' in samp_string:
                    pop = 'ancestor'
                    bc = 'N/A,N/A'
                else:
                    pop = evo_cond['population'][np.where(evo_cond['origin_track']==samp_string)[0][0]]
                    bc = barcode_dict[samp_string]

                chrom = variant.CHROM
                position = int(variant.end)
                mut_loc = f'{chrom}:{position}'
                gene = find_gene_at_position(genes, chrom, position - 1)  # Adjust for 0-based indexing
                if gene:
                    gene_name = gene['gene_id']
                else:
                    gene_name = 'NO_GENE'

                file.write(','.join([samp_string,mut_loc,gene_name,pop,bc,cov_dict[samp_string],str(variant.POS)])+'\n')

                # write sites with more than one different mutation to an additional other file
                if len(samp_list)>1:
                    multip_file.write(','.join([samp_string,mut_loc,gene_name,pop,bc,cov_dict[samp_string],str(variant.POS)])+'\n')
                snp_count += 1

                vcf_writer.write_record(variant)

print('number of variable sites after filtering: %s'%snp_count)

with open(filtered_file,'w') as file:
    # for variant in hard_filtered:
    #     chrom = variant.CHROM
    #     position = int(variant.end)
    #     mut_loc = f'{chrom}:{position}'
    #     file.write(mut_loc+',hard_filtered\n')

    for variant in ancestral_variants:
        chrom = variant.CHROM
        position = int(variant.end)
        mut_loc = f'{chrom}:{position}'
        file.write(mut_loc+',presumed ancestral\n')

    for var_tuplet in multiplet_filtered:
        variant,samp_list = var_tuplet
        chrom = variant.CHROM
        position = int(variant.end)
        mut_loc = f'{chrom}:{position}'
        file.write(f'{mut_loc},occurred in {len(samp_list)} samples including '+','.join(samp_list[:2])+'\n')

    for variant in insufficient_var_depth:
        chrom = variant.CHROM
        position = int(variant.end)
        mut_loc = f'{chrom}:{position}'
        file.write(f'{mut_loc},occurred in no samples with sufficient depth\n')
