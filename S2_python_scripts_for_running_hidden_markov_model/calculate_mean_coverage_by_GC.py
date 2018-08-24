#!/usr/bin/python

# This script calculates the mean coverage by GC window for a list of samples from Ag1000G

import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
import numpy as np
import pandas as pd
from re import *

# Write a function to normalised the coverage for each GC content bin
def normalise_coverage_by_GC(coverage, mean_coverage_by_GC, ploidy = 2):
	output = coverage.copy()
	# For each counts value in the output object, associate it with the mean coverage for its GC bin
	output['expcov'] = [mean_coverage_by_GC.loc[x] for x in output['GC']]
	# Now divide each counts value by the coverage associated with its GC bin, we multiply by 2 because
        # that was the ploidy of our mean coverage.
	output['Normalised_coverage'] = ploidy * output['Counts'] / output['expcov']
	# In very rare cases, the expected coverage ends up being 0 (because none of the few windows with a given
	# GC content have any coverage). In these cases, we manually set the coverage to be 0 (rather than NaN)
	which_zeros = (output['expcov'] == 0)
	output.loc[(which_zeros, 'Normalised_coverage')] = 0
	# Return the new dataframe 
	return output

# Check how many arguments you have and do anything necessary
if (len(argv) == 4):
	base_folder = argv[1]
	input_gc_filename = argv[2]
	accessibility_threshold = float(argv[3])
else:
	raise Exception("Fail. There should be three command line arguments (base_folder, input_gc_filename, accessibility_threshold)")

print 'Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n'
print 'Input arguments:\n\tbase_folder: ' + base_folder + '\n\tinput_gc_filename: ' + input_gc_filename + '\n\taccessibility_threshold: ' + str(accessibility_threshold) + '\n\n'
stdout.flush()

# Set the chromosomes and the autosomes
chroms = ['2L', '2R', '3L', '3R', 'X']
sex_chrom = 'X'

sample_manifest_file = file(base_folder + '/phase2.AR1/samples/samples.manifest', 'r')
sample_names = [x.rstrip('\n') for x in sample_manifest_file.readlines()]
sample_names = [sub('-', '_', x) for x in sample_names]
sample_manifest_file.close()

# Load up the GC composition file
gc_all = pd.read_csv(input_gc_filename, sep = '\t', header = None)
# Round the GC content to the nearest integer percentage
gc_all.iloc[:,2] = (gc_all.iloc[:,2]*100).astype(int)
gc_grouped = gc_all.groupby(0)

# Load up the accessibility data
acc_all = pd.read_csv(base_folder + '/phase2.AR1/accessibility/mean_accessibility.csv', sep = '\t')
acc_grouped = acc_all.groupby('Chrom')

# Get the distribution of unmasked GC bins
gc_mainchroms = gc_all.loc[gc_all[0].apply(lambda x:x in chroms), :]
gc_unmasked = gc_mainchroms.loc[np.array(acc_all['Mean_accessibility'] >= accessibility_threshold), :]
output_table = pd.DataFrame(gc_unmasked.groupby(2).size(), columns = ['bin_freq'])
output_table.index = output_table.index.rename('GC')

# Check that the GC content and accessibility have the same number of bins for each chromosome
for chrom in chroms:
	if gc_grouped.get_group(chrom).shape[0] != acc_grouped.get_group(chrom).shape[0]:
		raise Exception('Fail. GC content and accessibility on chromosome ' + chrom + ' should have the same number of bins.')

# Here is where we will store the variance data
output_variance = pd.DataFrame(0, columns = chroms + ['autosomes'], index = sample_names)

# Now for each sample combine the counts data for all the autosomes and calculate the mean
for this_sample in sample_names:
	print '\nAnalysing sample ' + this_sample + '.'
	stdout.flush()
	# create the filenames for each autosome
	for chrom in chroms:
		# Get the counts data
		this_filename = base_folder + '/CNV_v2/counting_output_v4_' + chrom + '/phase2/counts_for_mHMM_python2_v4_' + this_sample + '_' + chrom + '_nomapq_output.txt'
		print '\tLoading file ' + this_filename 
		stdout.flush()
		these_counts = pd.read_csv(this_filename, sep = '\t')
		this_gc = gc_grouped.get_group(chrom)
		this_acc = acc_grouped.get_group(chrom)
		# Check that the number of bins for the coverage is correct (we have already checked that accessibility
		# and GC have the same number of bins).
		if these_counts.shape[0] != this_gc.shape[0]:
			raise Exception('Fail. Coverage and GC_content on chromosome ' + chrom + ' should have the same number of bins.')
		# Add the GC content information to the counts table. 
		these_counts['GC'] = this_gc.iloc[:,2].values
		these_counts['Chrom'] = chrom
		# Get the accessibility-masked coverage
		these_masked_counts = these_counts.loc[np.array(this_acc['Mean_accessibility'] >= accessibility_threshold), :]
		# Add the data to the autosomal or sex chromsome record
		if chrom == sex_chrom:
			sexchrom_masked_counts = these_masked_counts
		else:
			if chrom == '2L':
				autosomal_masked_counts = these_masked_counts
			else:
				autosomal_masked_counts = pd.concat([autosomal_masked_counts, these_masked_counts])
			
	# Now we have combined the results from all of the autosomes, we can group the coverage by GC bin
	# and then compute the mean
	counts_by_GC = autosomal_masked_counts.groupby('GC')
	mean_counts_by_GC = counts_by_GC['Counts'].mean()
	# Add these to the global data
	output_table[this_sample] = mean_counts_by_GC

	# Now we want to calculate the variance in coverage for this sample, both for each chromosome
	# and overall
	print '\tCalculating normalised coverage and variance.'
	stdout.flush()
	masked_counts_by_autosome = autosomal_masked_counts.groupby('Chrom')
	for chrom in chroms:
		if chrom == sex_chrom:
			this_normcov = normalise_coverage_by_GC(sexchrom_masked_counts, mean_counts_by_GC)
		else:
			this_normcov = normalise_coverage_by_GC(masked_counts_by_autosome.get_group(chrom), mean_counts_by_GC)
		output_variance.loc[this_sample, chrom] = np.var(this_normcov['Normalised_coverage'])
	# Now calculate the variance across the autosomes
	output_variance.loc[this_sample, 'autosomes'] = np.var(normalise_coverage_by_GC(autosomal_masked_counts, mean_counts_by_GC)['Normalised_coverage'])
	
		
# Write the output to file
output_filename = 'mean_coverage_by_GC_masked_' + sub('\.', '', str(accessibility_threshold)) + '.csv'
print '\nSaving coverage output to file ' + output_filename 
stdout.flush()
output_table.to_csv(output_filename, sep = '\t')

output_variance_filename = 'coverage_variance_masked_' + sub('\.', '', str(accessibility_threshold)) + '.csv'
print '\nSaving variance output to file ' + output_variance_filename 
stdout.flush()
output_variance.to_csv(output_variance_filename, sep = '\t')

print '\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n'

