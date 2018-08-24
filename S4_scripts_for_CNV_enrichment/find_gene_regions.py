#!/usr/bin/python
scriptname = "find_gene_regions.py"

# This script parses a base features file, find the genes and outputs them to a csv file

import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
import pandas as pd
from re import *

# Run the script:

if __name__ == '__main__':
	
	# Check how many arguments you have and do anything necessary
	if (len(argv) == 2):
		input_filename = argv[1]
		output_filename = 'find_gene_regions_output.csv'
		region_key = 'gene'
	elif (len(argv) == 3):
		input_filename = argv[1]
		output_filename = argv[2]
		region_key = 'gene'
	elif (len(argv) == 4):
		input_filename = argv[1]
		output_filename = argv[2]
		region_key = argv[3]
	else:
		raise Exception("Fail. There should be one to three command line arguments (input_filename [, output_filename [, region_key]])")


	print 'Running ' + scriptname  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n'
	stdout.flush()

	# Load up the input file
	features_file = open(input_filename, 'r')

	# Create the table that will store the results
	output_table = pd.DataFrame(columns = ['Chrom', 'start', 'end'])

	# Go through the file and, when you find a gene, save it to the table
	for line in features_file:
		if line[0] == '#':
			continue
		info = line.split('\t')
		if info[2] == region_key:
			gene_name_string = info[8]
			gene_name = findall('(?<=ID=).+?(?=;)', gene_name_string)[0]
			output_table.loc[gene_name] = [info[0]] + info[3:5]
	
	# Write output to file
	output_table.to_csv(output_filename, sep = '\t')

	print '\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n'
	stdout.flush()

