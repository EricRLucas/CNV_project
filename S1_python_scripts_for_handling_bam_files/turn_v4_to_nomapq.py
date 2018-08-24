#!/usr/bin/python

# This script turns the output files from the counts_for_mHMM_python2_v4.py script into one where the counts 
# column contains the sum of all columns

import pandas as pd
from sys import argv
from re import *

input_filenames_list = argv[1]
input_filenames = [x.rstrip('\n') for x in file(input_filenames_list, 'r').readlines()]

for filename in input_filenames:
    this_counts_file = pd.read_csv(filename, sep = '\t')
    summed_counts =  this_counts_file['Counts mapq >= 10'] + this_counts_file['Counts mapq0'] + this_counts_file['Counts 0 < mapq < 10']
    this_output_file = pd.DataFrame(this_counts_file['Position'])
    this_output_file['Counts'] = summed_counts
    output_filename = sub('_output.txt', '_nomapq_output.txt', filename)
    this_output_file.to_csv(output_filename, sep = '\t', index = False)

