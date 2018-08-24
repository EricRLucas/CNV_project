#!/usr/bin/python
#find_target_genes.py

# This script goes through a list of transcripts and identifies genes of potential interest based on 
# certain keywords. 

import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
from re import *


# Set the format to identify the transcript name from within the transcript description
gene_name_format = '(?<=^>).+(?=-R)'

# Run the script:

if (len(argv) < 4):
	raise Exception("Fail. There should at least three command line argument (input_filename, output_filename_root, *search terms*)")

input_filename = argv[1]
output_filename_root = argv[2]
# Here are the terms for genes we are interested in
search_terms = argv[3:]

print 'Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n'
stdout.flush()

print 'Input arguments:\n\tinput_filename: ' + input_filename + '\n\toutput_filename_root: ' + output_filename_root + '\n\tsearch terms: ' + '; '.join(search_terms) + '\n\n'
stdout.flush()

# Load the input file
inputfile = file(input_filename, 'r')

# Prepare the output files for writing
outputfile1 = file(output_filename_root + '_output1.txt', 'w')
outputfile2 = file(output_filename_root + '_output2.txt', 'w')

# This object will store all the different genes we have found so far
genes_found_so_far = []

# Go through every line of the input file and identify genes of interest
for line in inputfile:
	for term in search_terms:
		search_term = '(?i)^>.+' + term
		if search(search_term, line):
			outputfile1.write(line)
			gene_name = findall(gene_name_format, line)[0]
			if gene_name not in genes_found_so_far:
				genes_found_so_far += [gene_name]
				outputfile2.write(gene_name + '\n')
			break

outputfile1.close()
outputfile2.close()

# And that's the end of the script
print '\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n'
stdout.flush()


