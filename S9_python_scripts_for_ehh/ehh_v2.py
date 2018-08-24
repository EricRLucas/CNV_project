#ehh_v2.py

# This script takes a vcf file of CNVs to create groups of core haplotypes and then uses haplotype data
# to calculate EHH for each of those core haplotypes.

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
	chrom = argv[1]
	core_left = argv[2]
	core_right = argv[3]
	h5_filename = argv[4]
	CNV_vcf_filename = argv[5]
	sample_info_filename = argv[6]
	ehh_output_filename = argv[7]
	hapshare_output_filename = argv[8]
else:
	raise Exception("Fail. There should be eight command line arguments (chrom, core_left, core_right, h5_filename, CNV_vcf_filename, sample_info_filename, ehh_output_filename, hapshare_output_filename)")

print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')
stdout.flush()

print('Input arguments:\n\tchrom: ' + chrom + '\n\tcore_left: ' + core_left + '\n\tcore_right: ' + core_right + '\n\th5_filename: ' + h5_filename + '\n\tsample_info_filename: ' + sample_info_filename + '\n\tCNV_vcf_filename: ' + CNV_vcf_filename + '\n\tehh_output_filename: ' + ehh_output_filename + '\n\thapshare_output_filename: ' + hapshare_output_filename + '\n\n')
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

def extract_haplotype_array(haparray, haplotype_pos, core_left, core_right, flank):
	# Get the haplotype on the right of the group.
	loc_right = haplotype_pos.locate_range(core_right, core_right + flank)
	haps_right_3d = haparray[loc_right, :1142, :]
	haps_right = allel.HaplotypeArray(haps_right_3d.reshape(haps_right_3d.shape[0], 2284))
	pos_right = haplotype_pos[loc_right]
	# Get the haplotype on the left of the group.
	loc_left = haplotype_pos.locate_range(core_left - flank, core_left)
	haps_left_3d = haparray[loc_left, :1142, :]
	haps_left = allel.HaplotypeArray(haps_left_3d.reshape(haps_left_3d.shape[0], 2284))
	pos_left = haplotype_pos[loc_left]
	
	return (pos_left, haps_left, pos_right, haps_right)

# Write a function to obtain the length of the shared haplotype between all pairwise comparisons in a group
# of haplotypes. 
def pairwise_shared_hapl(haparray, haplotype_pos, core_left, core_right, core_haps, flank1, flank2):
	# Create two sets of data, one with the small flank range and one with the large flank range
	pos_left, haps_left, pos_right, haps_right = extract_haplotype_array(haparray, haplotype_pos, core_left, core_right, flank1)
	pos_left_bigflank, haps_left_bigflank, pos_right_bigflank, haps_right_bigflank = extract_haplotype_array(haparray, haplotype_pos, core_left, core_right, flank2)
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
				this_hap_pair_left = haps_left.take([index1,index2], axis = 1) 
				this_hap_pair_right = haps_right.take([index1,index2], axis = 1) 
				# Calculate the shared haplotype using the small flank value. If this is too small, (ie: if the 
				# two samples are identical at this range), try again with the large one.
				shared_hap_length_left = shared_hapl(this_hap_pair_left, pos_left, True)
				if shared_hap_length_left == pos_left[-1] - pos_left[0] + 1:
					this_hap_pair_left = haps_left_bigflank.take([index1,index2], axis = 1) 
					shared_hap_length_left = shared_hapl(this_hap_pair_left, pos_left_bigflank, True)
					if shared_hap_length_left == pos_left_bigflank[-1] - pos_left_bigflank[0] + 1:
						print('All SNPs were identical in this haplotype pair')
				shared_hap_length_right = shared_hapl(this_hap_pair_right, pos_right, False)
				if shared_hap_length_right == pos_right[-1] - pos_right[0] + 1:
					this_hap_pair_right = haps_right_bigflank.take([index1,index2], axis = 1) 
					shared_hap_length_right = shared_hapl(this_hap_pair_right, pos_right_bigflank, False)
					if shared_hap_length_right == pos_right_bigflank[-1] - pos_right_bigflank[0] + 1:
						print('All SNPs were identical in this haplotype pair')
				all_pairwise_lengths_left += [shared_hap_length_left]
				all_pairwise_lengths_right += [shared_hap_length_right]
		output_dict[this_name] = (all_pairwise_lengths_left, all_pairwise_lengths_right)
	return(output_dict)
	
# Write a function to get the EHH decay around a core haplotype
def get_ehh_decay(haparray, haplotype_pos, core_left, core_right, core_haps, flank):
	# Get the haplotype arrays
	pos_left, haps_left, pos_right, haps_right = extract_haplotype_array(haparray, haplotype_pos, core_left, core_right, flank)
	# Get the positions of interest
	all_pos = np.concatenate((pos_left, pos_right), axis = 0)
	# Calculate the EHH
	all_ehh = pd.DataFrame(columns = all_pos)
	for (l, s) in core_haps:
		haps_left_core = haps_left.take(sorted(s), axis=1)
		haps_right_core = haps_right.take(sorted(s), axis=1)
		ehh_decay_left = allel.ehh_decay(haps_left_core[::-1])
		ehh_decay_right = allel.ehh_decay(haps_right_core)
		all_ehh.loc[l,:] = np.concatenate((ehh_decay_left[::-1], ehh_decay_right), axis = 0)
	
	return all_ehh.transpose()

# Load the callset
print('Loading CNV vcf file.')
stdout.flush()
callset = allel.read_vcf(CNV_vcf_filename)
# Get the metadata
sample_info = pd.read_csv(sample_info_filename, sep = '\t', index_col = 0)
# Get the population information and repeat each entry twice to have a record of haplotype populations
haplotype_pop = []
for p in sample_info['population']:
	haplotype_pop.extend([p, p])
# Get the CNV_calls. These are initially produced as a 3-dimensional array and are then flattened to 2
# dimensions. 
CNV_calls_3d = callset['calldata/GT']
CNV_calls = CNV_calls_3d.reshape(CNV_calls_3d.shape[0], 2284)
# For each haplotype index, we want to know which CNV it carries, if any. We do this by 
# pasting together all of the CNV IDs for that haplotype
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
print('Loading haplotypes.')
stdout.flush()
hapsfile = h5py.File(h5_filename)
allhaps = hapsfile[chrom]['calldata']['genotype']
pos = allel.SortedIndex(hapsfile[chrom]['variants']['POS'])

# Get the ehh decay in a table
print('Calculating EHH in each group.')
stdout.flush()
output_table = get_ehh_decay(allhaps, pos, int(core_left), int(core_right), CNV_core_tuple_splitwt, 100000)

# Now write the output table to file
print('Writing EHH output to file.')
stdout.flush()
output_table.to_csv(ehh_output_filename, sep = '\t')

print('Calculating length of shared haplotypes in each group.')
stdout.flush()
# For each core haplotype with more than two samples, calculate the shared haplotype lengths between 
# all pairs of haplotypes.
hapshare = pairwise_shared_hapl(allhaps, pos, int(core_left), int(core_right), CNV_core_tuple_splitwt, 100000, 4000000)

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

