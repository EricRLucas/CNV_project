#!/usr/bin/python

# This script is modified from the one found here: http://nbviewer.jupyter.org/github/alimanfoo/hmmcnv/blob/master/tutorial.ipynb
# It runs an mHMM on coverage counts obtained from a bamfile, using normalisation based on the mean coverage per 
# GC bin over the whole genome of the individual.

import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
from hmmlearn.hmm import GaussianHMM
import numpy as np
import pandas as pd
from re import *

def fit_hmm(depth_normed,  # normalised coverage array 
            transition_probability,  # probability of state transition
            variance,  # variance per copy 
            variance_fixed,  # variance for the zero copy number state 
            max_copy_number=12,  # maximum copy number to consider in the model 
            n_iter=0,  # number of iterations to perform when fitting the model
            params='st',  # parameters that can be changed through fitting 
            init_params=''  # parameters that are initialised from the data
           ):
    
	# convenience variable
	min_copy_number = 0  # minimum copy number to consider in the model
	n_states = max_copy_number - min_copy_number + 1
	
	# construct the transition matrix
	transmat = np.zeros((n_states, n_states))
	transmat[:] = transition_probability
	transmat[np.diag_indices(n_states)] = 1-((n_states-1)*transition_probability)
	
	# construct means and covariance
	means_list = range(n_states)
	means = np.array([[n] for n in means_list])
	covars = np.array([[variance*n + variance_fixed] for n in means_list])
	
	# setup HMM 
	model = GaussianHMM(n_states, 
	                    covariance_type='diag', 
	                    n_iter=n_iter, 
	                    params=params,
	                    init_params=init_params)
	model.means_ = means
	model.covars_ = covars
	model.transmat_ = transmat
	
	# fit HMM
	obs = np.column_stack([depth_normed])
	model.fit(obs)
	
	# predict hidden states
	h = model.predict(obs)
	
	return h

# Write a function to normalised the coverage for each GC content bin
def normalise_coverage_by_GC(coverage, mean_coverage_by_GC, ploidy = 2):
	output = coverage.copy()
	# For each counts value in the output object, associate it with the mean coverage for its GC bin
	output['expcov'] = [mean_coverage_by_GC.loc[x] for x in output['GC']]
	# Now divide each counts value by the coverage associated with its GC bin, we multiply by 2 because
	# that was the ploidy of the expected coverage.  
	output['Normalised_coverage'] = ploidy * output['Counts'] / output['expcov']
	# Return the new dataframe 
	return output


# Check how many arguments you have and do anything necessary
if (len(argv) == 7):
	input_counts_filename_list = argv[1]
	input_gc_filename = argv[2]
	mean_coverage_by_gc_filename = argv[3]
	variance_by_sample_filename = argv[4]
	mapq_proportions_filename = argv[5]
	max_mapq = argv[6]
	output_filename_root = 'mHMM_output.txt'
elif (len(argv) == 8):
	input_counts_filename_list = argv[1]
	input_gc_filename = argv[2]
	mean_coverage_by_gc_filename = argv[3]
	variance_by_sample_filename = argv[4]
	mapq_proportions_filename = argv[5]
	max_mapq = argv[6]
	output_filename_root = argv[7]
else:
	raise Exception("Fail. There should be six or seven command line arguments (input_counts_filename_list, input_gc_filename, mean_coverage_by_gc_filename, variance_by_sample_filename, mapq_proportion_filename, max_mapq, [, output_filename_root])")


print 'Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n'
print 'Input arguments:\n\tinput_counts_filename_list: ' + input_counts_filename_list + '\n\tinput_gc_filename: ' + input_gc_filename + '\n\tmean_coverage_by_gc_filename: ' + mean_coverage_by_gc_filename + '\n\tvariance_by_sample_filename: ' + variance_by_sample_filename + '\n\tmapq_proportions_filename: ' + mapq_proportions_filename + '\n\tmax_mapq: ' + max_mapq + '\n\toutput_filename_root: ' + output_filename_root + '\n\n'
stdout.flush()

max_mapq = float(max_mapq)
# Load up the file containing the mapq information and calculate the proportion of mapq0 reads in each bin
mapq_prop = pd.read_csv(mapq_proportions_filename, sep = '\t') 
mapq_prop['mapq0_prop'] = mapq_prop['Count mapq = 0'] / (mapq_prop['Count mapq = 0'] + mapq_prop['Count mapq > 0'])
# Create a boolean for whether the mapq0 proportion is low enough
mapq_prop['mapq0_acceptable'] = mapq_prop['mapq0_prop'] <= max_mapq
# Split the data by chromosome
mapq_prop_bychrom = mapq_prop.groupby('Chrom')


# Load up the list of counts filenames
list_of_counts_filenames = file(input_counts_filename_list, 'r').readlines()
# Remove the trailing newline character from each filename
list_of_counts_filenames = [x.rstrip('\n') for x in list_of_counts_filenames]

# Get the GC content
gc_all = pd.read_csv(input_gc_filename, sep = '\t', header = None)
# Round the GC content to the nearest integer percentage
gc_all.iloc[:,2] = (gc_all.iloc[:,2]*100).astype(int)
# Load up the table of mean coverage per GC
mean_cov_by_gc_all_bins = pd.read_csv(mean_coverage_by_gc_filename, sep = '\t', index_col = 'GC')
# Remove bins that don't have at least n = 100 in the accessible genome
mean_cov_by_gc = mean_cov_by_gc_all_bins.loc[mean_cov_by_gc_all_bins['bin_freq'] >= 100, :]

# Load up the table of mean variance per sample
var_by_sample = pd.read_csv(variance_by_sample_filename, sep = '\t', index_col = 0)

# Go through each of the files
for input_counts_filename in list_of_counts_filenames:
	print '\nLoading file ' + input_counts_filename
	stdout.flush()
	# Load up the counts data
	input_coverage = pd.read_csv(input_counts_filename, sep = '\t')

	# Determine the sample name and chromsome from the counts file
	sample_name = findall('[A-Z]{2}\d{4}_C[^_]?', input_counts_filename)[0] 
	# Because look behind requires fixed width parameters, we need to look for the chromosome in two stages
	prechrom = findall('(?<=[A-Z]{2}\d{4}_C)[^_]?_[^_]+', input_counts_filename)[0] 
	chrom = findall('[^_]+$', prechrom)[0]

	# Get the GC content for this chromosome
	gc_thischrom = gc_all.groupby(0).get_group(chrom)
	# Get the mapq0 information for this chromosome
	mapq_prop_thischrom = mapq_prop_bychrom.get_group(chrom)
	# Check that the length of the coverage and GC_content dataframes are the same
	if input_coverage.shape[0] != gc_thischrom.shape[0]:
		raise Exception('Fail. Coverage and GC_content should have the same number of elements')
	# Check that the length of the coverage and mapq0 dataframes are the same
	if input_coverage.shape[0] != mapq_prop_thischrom.shape[0]:
		raise Exception('Fail. Coverage and GC_content should have the same number of elements')

	# We add the GC content data to the counts table 
	input_coverage['GC'] = gc_thischrom.iloc[:,2].values

	# Get the mean coverage by GC bin for this sample
	this_sample_mean_cov = mean_cov_by_gc[sample_name]

	# Mask out the positions with low mapq0.
	input_coverage_mapq_masked = input_coverage.loc[np.array(mapq_prop_thischrom['mapq0_acceptable']),:]
	print str(np.sum(-mapq_prop_thischrom['mapq0_acceptable'])) + ' out of ' + str(input_coverage.shape[0]) + ' were masked for having more than ' + str(max_mapq * 100) + '% mapq0 reads.'
	# Mask out the positions that don't have a GC bin with mean coverage data. 
	input_coverage_mapqandgc_masked = input_coverage_mapq_masked.loc[input_coverage_mapq_masked['GC'].isin(this_sample_mean_cov.index),:]
	print 'Of these, ' + str(np.sum(-input_coverage_mapq_masked['GC'].isin(this_sample_mean_cov.index))) + ' were masked for having a GC bin with too little coverage information.'
	
	# Normalised coverage according to the GC content
	print '\tNormalising coverage'
	stdout.flush()
	main_table = normalise_coverage_by_GC(input_coverage_mapqandgc_masked, this_sample_mean_cov)

	# Get the autosomal variance in coverage for this sample
	this_mean_variance = var_by_sample.loc[(sample_name, 'autosomes')]

	# Fit the CNV model. The variance we have calculated is for the diploid state, but the model expects
	# the variance increment unit. So we feed it the variance / 2. 
	trans_prob = 0.00001
	print '\tFitting CNV model with trans_prob = ' + str(trans_prob) + '.'
	stdout.flush()
	cnv = fit_hmm(main_table['Normalised_coverage'].values,  
				  transition_probability = trans_prob,
				  variance = this_mean_variance/2,  
				  variance_fixed = 0.001, 
				  max_copy_number=12)
	# Add the CNV values to the table
	main_table['CNV'] = cnv

	# Build a unique output filename for this sample
	output_filename = sub('\.[^.]+$', '', sub('.*/', '', input_counts_filename)) + '_' + output_filename_root + '.csv'
	# Write that table to file
	print '\tSaving output to file ' + output_filename
	stdout.flush()
	main_table.to_csv(output_filename, sep = '\t')

print '\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n'
stdout.flush()

