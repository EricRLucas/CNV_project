#ehh_X_v2.py

# This script takes a vcf file of CNVs to create groups of core haplotypes and then uses haplotype data
# to calculate EHH for each of those core haplotypes. It specifically works with data from the X chromosome,
# and therefore needs genotype as well as haplotype information (because the haplotype data does not include
# data from males). 

import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
from re import *
import allel
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict as dict 
from re import *

if (len(argv) == 9):
	core_left = argv[1]
	core_right = argv[2]
	haplotype_filename = argv[3]
	genotype_filename = argv[4]
	CNV_vcf_filename = argv[5]
	sample_info_filename = argv[6]
	ehh_output_filename = argv[7]
	hapshare_output_filename = argv[8]
else:
	raise Exception("Fail. There should be seven command line arguments (core_left, core_right, haplotype_filename, genotype_filename, CNV_vcf_filename, sample_info_filename, ehh_output_filename, hapshare_output_filename)")

print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')
stdout.flush()

print('Input arguments:\n\tcore_left: ' + core_left + '\n\tcore_right: ' + core_right + '\n\thaplotype_filename: ' + haplotype_filename + '\n\tgenotype_filename: ' + genotype_filename + '\n\tsample_info_filename: ' + sample_info_filename + '\n\tCNV_vcf_filename: ' + CNV_vcf_filename + '\n\tehh_output_filename: ' + ehh_output_filename + '\n\thapshare_output_filename: ' + hapshare_output_filename + '\n\n')
stdout.flush()

# Create a function to calculate the length of the shared haplotype between a pair of samples. The hap_pair 
# object is an array with just two columns, where each column represents a haplotype. 
def shared_hapl(hap_pair, positions, leftside = False):
	# Check that the number of columns in hap_pair is 2
	if hap_pair.shape[1] != 2:
		raise Exception('hap_pair should have exactly two columns.')
	# If this is to the left of the core, invert the haplotype array
	if leftside:
		hap_pair = hap_pair[::-1]
	# Create a list of booleans comparing the two values in each row
	same_hap = np.where(hap_pair[:,0] != hap_pair[:,1])[0]
	# If same_hap does not exist, then the haplotypes are identical. This implies that we should use a larger haplotype
	# range as we are only getting a minimum shared haplotype estimate
	if len(same_hap) == 0:
		shared_distance = positions[-1] - positions [0] + 1
	else :
		# The first different position is where the haplotype breaks down. The position before this, minus the start position, 
		# is the length of the shared haplotype. If the first position is already different, we give it a shared haplotype of
		# 0. 
		last_shared_index = min(same_hap) - 1
		if last_shared_index == -1:
			shared_distance = 0
		else:
			if leftside:
				shared_distance = positions[-1] - positions[-last_shared_index - 1] + 1
			else:
				shared_distance = positions[last_shared_index] - positions[0] + 1
	return(shared_distance)

def extract_haplotype_array(haparray, genarray, haplotype_pos, genotype_pos, male_indices, core_left, core_right, flank):
	# Get the haplotype on the right of the group.
	haps_loc_right = haplotype_pos.locate_range(core_right, core_right + flank)
	haps_right_3d = haparray[haps_loc_right, :1058, :]
	haps_right = allel.HaplotypeArray(haps_right_3d.reshape(haps_right_3d.shape[0], 1058*2))
	pos_right = haplotype_pos[haps_loc_right]
	# Get the genotypes on the right of the group. Some genotypes are not present in the haplotype table so
	# we need to remove them. It takes a long time to do this directly, and is much quicker to first take a
	# slice for the correct range, and then filter out the loci not found in the haplotypes data. 
	gen_loc_right = genotype_pos.locate_range(core_right, core_right + flank)
	gen_right_all_3d = genarray[gen_loc_right, :, :]
	gen_pos_right = genotype_pos[gen_loc_right]
	gen_loc_right_inhap = gen_pos_right.locate_keys(pos_right)
	gen_right_all_inhap_3d = gen_right_all_3d[gen_loc_right_inhap, :, :]
	gen_right = allel.HaplotypeArray(gen_right_all_inhap_3d[:, male_indices, 0])
	# Combine the haplotypes and genotypes
	genhaps_right_all = allel.HaplotypeArray(np.concatenate((haps_right, gen_right), 1))
	# Get the haplotype on the left of the group.
	haps_loc_left = haplotype_pos.locate_range(core_left - flank, core_left)
	haps_left_3d = haparray[haps_loc_left, :1058, :]
	haps_left = allel.HaplotypeArray(haps_left_3d.reshape(haps_left_3d.shape[0], 1058*2))
	pos_left = haplotype_pos[haps_loc_left]
	# Get the genotypes on the left of the group
	gen_loc_left = genotype_pos.locate_range(core_left - flank, core_left)
	gen_left_all_3d = genarray[gen_loc_left, :, :]
	gen_pos_left = genotype_pos[gen_loc_left]
	gen_loc_left_inhap = gen_pos_left.locate_keys(pos_left)
	gen_left_all_inhap_3d = gen_left_all_3d[gen_loc_left_inhap, :, :]
	gen_left = allel.HaplotypeArray(gen_left_all_inhap_3d[:, male_indices, 0])
	# Combine the haplotypes and genotypes
	genhaps_left_all = allel.HaplotypeArray(np.concatenate((haps_left, gen_left), 1))
	# We need to exclude loci where there are missing values (which occurs in the genotype table)
	no_missing_values_right = np.apply_along_axis(lambda x: len(np.where(x == -1)[0]) == 0, 1, genhaps_right_all)
	genhaps_right = genhaps_right_all[no_missing_values_right, :]
	no_missing_values_left = np.apply_along_axis(lambda x: len(np.where(x == -1)[0]) == 0, 1, genhaps_left_all)
	genhaps_left = genhaps_left_all[no_missing_values_left, :]
	
	return (pos_left[no_missing_values_left], genhaps_left, pos_right[no_missing_values_right], genhaps_right)

# Write a function to obtain the length of the shared haplotype between all pairwise comparisons in a group
# of haplotypes. 
def pairwise_shared_hapl(haparray, genarray, haplotype_pos, genotype_pos, male_indices, core_left, core_right, core_haps, flank1, flank2):
	# Create two sets of data, one with the small flank range and one with the large flank range
	pos_left, genhaps_left, pos_right, genhaps_right = extract_haplotype_array(haparray, genarray, haplotype_pos, genotype_pos, male_indices, core_left, core_right, flank1)
	pos_left_bigflank, genhaps_left_bigflank, pos_right_bigflank, genhaps_right_bigflank = extract_haplotype_array(haparray, genarray, haplotype_pos, genotype_pos, male_indices, core_left, core_right, flank2)
	# Create a dictionary to store the output
	output_dict = dict()
	for this_name, this_core in core_haps:
		print(this_name)
		stdout.flush()
		if len(this_core) < 2:
			output_dict[this_name] = ([None], [None])
			continue
		# We want to get all of the pairwise combinations of haplotypes in this group
		all_pairwise_lengths_left = []
		all_pairwise_lengths_right = []
		for i in range(len(this_core)):
			for j in range(i+1,len(this_core)):
				index1 = this_core[i]
				index2 = this_core[j]
				this_hap_pair_left = genhaps_left.take([index1,index2], axis = 1) 
				this_hap_pair_right = genhaps_right.take([index1,index2], axis = 1) 
				# Calculate the shared haplotype using the small flank value. If this is too small, (ie: if the 
				# two samples are identical at this range), try again with the large one.
				shared_hap_length_left = shared_hapl(this_hap_pair_left, pos_left, True)
				if shared_hap_length_left == pos_left[-1] - pos_left[0] + 1:
					this_hap_pair_left = genhaps_left_bigflank.take([index1,index2], axis = 1) 
					shared_hap_length_left = shared_hapl(this_hap_pair_left, pos_left_bigflank, True)
					if shared_hap_length_left == pos_left_bigflank[-1] - pos_left_bigflank[0] + 1:
						print('All SNPs were identical in this haplotype pair')
				shared_hap_length_right = shared_hapl(this_hap_pair_right, pos_right, False)
				if shared_hap_length_right == pos_right[-1] - pos_right[0] + 1:
					this_hap_pair_right = genhaps_right_bigflank.take([index1,index2], axis = 1) 
					shared_hap_length_right = shared_hapl(this_hap_pair_right, pos_right_bigflank, False)
					if shared_hap_length_right == pos_right_bigflank[-1] - pos_right_bigflank[0] + 1:
						print('All SNPs were identical in this haplotype pair')
				all_pairwise_lengths_left += [shared_hap_length_left]
				all_pairwise_lengths_right += [shared_hap_length_right]
		output_dict[this_name] = (all_pairwise_lengths_left, all_pairwise_lengths_right)
	return(output_dict)

# Write a function to get the EHH decay around a core haplotype
def get_ehh_decay_X(haparray, genarray, haplotype_pos, genotype_pos, male_indices, core_left, core_right, core_haps, flank):
	# Get the haplotype arrays
	pos_left, genhaps_left, pos_right, genhaps_right = extract_haplotype_array(haparray, genarray, haplotype_pos, genotype_pos, male_indices, core_left, core_right, flank)
	# Get the positions of interest
	all_pos = np.concatenate((pos_left, pos_right), axis = 0)
	# Calculate the EHH
	all_ehh = pd.DataFrame(columns = all_pos)
	for (l, s) in core_haps:
		genhaps_left_core = genhaps_left.take(sorted(s), axis=1)
		genhaps_right_core = genhaps_right.take(sorted(s), axis=1)
		ehh_decay_left = allel.ehh_decay(genhaps_left_core[::-1])
		ehh_decay_right = allel.ehh_decay(genhaps_right_core)
		all_ehh.loc[l,:] = np.concatenate((ehh_decay_left[::-1], ehh_decay_right), axis = 0)
	
	return all_ehh.transpose()

chrom = 'X'

# Load the callset
print('Loading CNV vcf file.')
stdout.flush()
callset = allel.read_vcf(CNV_vcf_filename)
# Get the metadata
sample_info = pd.read_csv(sample_info_filename, sep = '\t', index_col = 0)
# Get the population information. For females,  repeat each entry twice to have a record of haplotype 
# populations. Do the females first and then the males to match the order of samples in the haplotype 
# and CNV data. Also record the sample names as you go along for a sanity check that the final order of 
# samples matches between the population data and the rest.
haplotype_fem_pop = []
haplotype_male_pop = []
fem_sample_names = []
male_sample_names = []
for s in sample_info.index:
	p = sample_info.loc[s, 'population']
	if sample_info.loc[s, 'sex'] == 'F':
		haplotype_fem_pop.extend([p, p])
		fem_sample_names.extend([s])
	else:
		haplotype_male_pop.extend([p])
		male_sample_names.extend([s])
# Join the outputs together
haplotype_pop = np.array(haplotype_fem_pop + haplotype_male_pop)
meta_samples = np.array(fem_sample_names + male_sample_names)
# As a sanity check, make sure that the samples are now the same as the CNV samples
CNV_call_samples = callset['samples']
if not np.array_equal(meta_samples, CNV_call_samples):
	raise Exception('Order of samples is not the same in the metadata and the ordered CNV calls.')

# Get the CNV_calls. These are initially produced as a 3-dimensional array and are then flattened to 2
# dimensions. 
CNV_calls_3d = callset['calldata/GT']
# Split these by males and females
CNV_calls_fem_3d = CNV_calls_3d[:,:1058,:]
CNV_calls_male_3d = CNV_calls_3d[:,1058:,:]
# For females, we flatten the 3d array into a 2d one by putting the pairs of haplotypes next to each other.
CNV_calls_fem = CNV_calls_fem_3d.reshape(CNV_calls_3d.shape[0], 1058*2)
# For males, we just take the first column, since the second column is just a repeat of the first. 
CNV_calls_male = CNV_calls_male_3d[:,:,0]
# Join the tables together
CNV_calls = np.concatenate((CNV_calls_fem, CNV_calls_male), 1)

# For each haplotype index, we want to know which CNV it carries, if any. We do this by pasting together 
# all of the CNV IDs for that haplotype
def paste_CNV_names(x, names):
	hapname = '_'.join(names[x.astype('bool')])
	if hapname == '':
		return 'wt'
	else:
		return hapname
CNV_cores = []
for i in range(CNV_calls.shape[1]):
	CNV_cores.append(paste_CNV_names(CNV_calls[:,i], callset['variants/ID']))
CNV_cores = np.array(CNV_cores)

# Get the CNVs we are primarily interested in
CNV_core_tuple = list()
allcores = np.unique(CNV_cores)
for cnv in allcores:
	CNV_core_tuple.append((cnv, np.where(CNV_cores == cnv)[0]))

# We want to split the wt core by population
wt_core = [x for x in CNV_core_tuple if x[0] == 'wt'][0][1]
wt_core_pops = np.array(haplotype_pop)[wt_core]
wt_core_tuple = list()
for p in np.unique(wt_core_pops):
	wt_core_tuple.append(('wt_' + p, wt_core[np.where(wt_core_pops == p)[0]]))
CNV_core_tuple_splitwt = wt_core_tuple + [x for x in CNV_core_tuple if x[0] != 'wt']

# Get the haplotypes
print('Loading female haplotypes.')
stdout.flush()
hapsfile = h5py.File(haplotype_filename)
allhaps = hapsfile[chrom]['calldata']['genotype']
pos = allel.SortedIndex(hapsfile[chrom]['variants']['POS'])

# Get the genotypes
print('Loading genotypes.')
stdout.flush()
gensfile = h5py.File(genotype_filename)
allgens = gensfile[chrom]['calldata']['genotype']
# To keep only the males, work out the index that will be needed to get them, and do a sanity check that 
# this is selecting the correct data. We don't reduce the table to just the males yet because we need to 
# first narrow down the range of loci involved, otherwise we get an enormous table. 
males_index = np.where(sample_info['sex'] == 'M')[0]
male_samples = [str(x, 'utf-8') for x in gensfile[chrom]['samples'].value[males_index]]
if not male_samples == male_sample_names:
	raise Exception('Failed to correctly retrieve the male samples from the genotype data.')
genpos = allel.SortedIndex(gensfile[chrom]['variants']['POS'])

# Get the ehh decay in a table
print('Calculating EHH in each group.')
stdout.flush()
output_table = get_ehh_decay_X(allhaps, allgens, pos, genpos, males_index, int(core_left), int(core_right), CNV_core_tuple_splitwt, 100000)

# Now write the output table to file
print('Writing EHH output to file.')
stdout.flush()
output_table.to_csv(ehh_output_filename, sep = '\t')

print('Calculating length of shared haplotypes in each group.')
stdout.flush()
# For each core haplotype with more than two samples, calculate the shared haplotype lengths between 
# all pairs of haplotypes.
hapshare = pairwise_shared_hapl(allhaps, allgens, pos, genpos, males_index, int(core_left), int(core_right), CNV_core_tuple_splitwt, 100000, 4000000)

# Now write this to the output file
print('Writing pairwise shared haplotype output to file.')
stdout.flush()
outputfile = open(hapshare_output_filename, 'w')
for this_name, (hap_left, hap_right) in hapshare.items():
	outputfile.write(this_name + '_left:\t' + ' '.join([str(x) for x in hap_left]) + '\n')
	outputfile.write(this_name + '_right:\t' + ' '.join([str(x) for x in hap_right]) + '\n')
outputfile.close()

print('\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n')
stdout.flush()

