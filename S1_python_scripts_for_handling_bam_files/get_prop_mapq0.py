#!/usr/bin/python3
scriptname = "get_prop_mapq0.py"

# This script takes the output to counts_for_mHMM_python2_v4.py for all individuals and calculates the 
# proportion of reads that have mapq=0

import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
from re import *

# Run the script:

if __name__ == '__main__':
	
	if (len(argv) == 3):
		list_of_filenames = argv[1]
		chromosome = argv[2]
		window_size = 300
		output_filename = 'get_prop_mapq0.txt'
	elif (len(argv) == 4):
		list_of_filenames = argv[1]
		chromosome = argv[2]
		window_size = argv[3]
		output_filename = 'get_prop_mapq0.txt'
	elif (len(argv) == 5):
		list_of_filenames = argv[1]
		chromosome = argv[2]
		window_size = argv[3]
		output_filename = argv[4]
	else:
		raise Exception("Fail. There should be between two and four command line arguments (list_of_filenames, chromosome [, window_size [, output_filename]])")


	print('Running ' + scriptname  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')
	stdout.flush()
	
	# Get file objects for all of the files
	files_to_load = [x.rstrip('\n') for x in open(list_of_filenames, 'r').readlines()]
	file_objects = [open(x, 'r') for x in files_to_load]

	# For each file, skip the header line
	for this_file in file_objects:
		first_line = next(this_file)
		if first_line != 'Position\tCounts mapq >= 10\tCounts 0 < mapq < 10\tCounts mapq0\n':
			raise Exception("Fail. The header row doesn't look at expected.")

	# Go through all files one by one and calculate the proportion of mapq=0 reads, checking each time that the 
	# positions match
	expected_position = 0
	output_file = open(output_filename, 'w')
	output_file.write('Position\tCount mapq > 0\tCount mapq = 0\n')
	should_break = False
	while(1):
		total_mapqpositive = 0
		total_mapq0 = 0
		for this_file in file_objects:
			try:
				this_line = next(this_file).rstrip('\n')
			except StopIteration:
				should_break = True
				break
			this_info = this_line.split('\t')
			if int(this_info[0]) != expected_position:
				raise Exception("Unexpected position in file " + this_file.name + ". Expected " + str(expected_position) + " but observed " + str(this_info[0]) + ".")
			total_mapq0 += int(this_info[3])
			total_mapqpositive += (int(this_info[1]) + int(this_info[2]))

		if should_break:
			break
		output_file.write(str(expected_position) + '\t' + str(total_mapqpositive) + '\t' + str(total_mapq0) + '\n')
		expected_position += 300
	
	output_file.close()

	print('\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n')
	stdout.flush()

