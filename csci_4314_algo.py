########################################################################
# Final Project: CSCI-4314
    # Date: April 2017
    
    # Purpose: Benchmark sequence aligner against Bowtie
    # Create a local aligner to benchmark against Bowtie
    # Compares a list of small sequences to one large genome and returns all alignments
    # Finds local regions with the highest levels of conservation/similarity

    # Psuedo-Code: Updated version of the Smith Waterman algorithm
    
    # Format: -SF <sequence_filename> -S <sequence> -AF <align_filename>
	# Sequence file: File that contains all large genomes that will be aligned against
		# Sequence: The specific sequence within the file that will be aligned against
	# Alignment File: File that contains all the smaller sequences to be aligned
########################################################################

import argparse # parse given flag options
import re
import operator # allows for itemgetter for max value in a dictionary

match_score = 3
gap_score = -1

########################################################################
## SETTING UP THE DICTIONARIES FROM THE GIVEN FILES

def readingFileDict(filename):
	# reads in the file and returns a dictionay with headers and sequences: {header:sequence}
	fullList = []
	seqlist = []
	header_list = []
	nucleo_list = []

	append = fullList.append # avoid re-using append (to improve running time)
	with open(filename, "r") as given_file:
		seq = ''
		for line in given_file:
			line = line.rstrip('\n').replace(" ", "@").replace("\t", "@").replace("\r", "@") # replace spaces with known character and replace tabs
			if line.startswith('>'):
				line = line + ' '
				line = line.replace('>', ' >')
			append(line)
		seq = ''.join(fullList).upper()
		seqList = seq.split()
		# Pulls out the sequences and genomes
		# By removing any extranous puncutation with predicatble characters to be spliced

		# Returns the header sequence name
		append = header_list.append # avoid re-using append (to improve running time)
		for element in seqList:
			if ">" in element:
				element = element.replace("@", "")
				element = element.strip(">") # assumes all header/sequences starts with >
				append(element) # returns a list of  headers ['chrI', 'chrII', etc...]

		nucleo_list =  [x for x in seqList if '>' not in x] # returns a list of sequence ['ATC', 'TGGC', etc..]
		nucleo_list = map(str.upper, nucleo_list) # convert all sequences to upper case for consitency

		for i in range(len(nucleo_list)): # remove extra @ left behind by replace and replace n with A (n = any ATCG)
			if 'N' in nucleo_list[i]:
				new_value = nucleo_list[i].replace("N", "A")
				nucleo_list[i] = new_value
			if '@' in nucleo_list[i]:
				new_value = nucleo_list[i].replace("@", "")
				nucleo_list[i] = new_value
		seq_dict = seqDictPairs(header_list, nucleo_list) # tuples of a pair's list and a dictionary {seq:gen}
	return seq_dict

def seqDictPairs(header_list, sequence_list):
	# creates a dictionary between the sequence (header) and the associated genome {seq:genome} dictionary
	# takes in the header list and the sequence list to combine into a single diectionary that can be searched
	# setting up sequence genome dictionary
	seq_gen_dict = {}
	seq_gen_dict = zip(header_list, sequence_list) # combine the two lists
	seq_gen_dict = dict(seq_gen_dict) # create new dictionary from the lists
	return seq_gen_dict

########################################################################
## CHECKS IF IT EXISTS EXACTLY (for optimizing)
def checkExists(align_seq, genome_seq):
	# return true/false if the expected sequence exists exactly (no gaps required)
	# returns (true/false, # of times it appears (default=0), index if it only appears once (default=None)
	pass

########################################################################
## CREATE MATRICES FOR LOCAL ALIGNMENT's DYNAMIC PROGRAMMING

def orderedPairs(align_dict, genome_seq):
	# creates pairs for each aligners and the genome
	# SEQ_0 and GEN_0, SEQ_1 and GEN_0
	total_pairs = []
	for aligner_sequence in align_dict:
		total_pairs.append([aligner_sequence, genome_seq])
		'''example: [['SEQ_1', 'GEN_1'], ['SEQ_0', 'GEN_1'], ['SEQ_3', 'GEN_1'], ['SEQ_2', 'GEN_1']]'''
	return total_pairs

def zeroMatrix(width_given, height_given):
	# sets up an matrix with zeros for the size of the sequences given, with gaps
	width, height = width_given, height_given
	matrix = [[0 for x in range(width)] for y in range(height)]
	return matrix

def DPlocalMatrix(sequence_large_genome, sequence_small_align, zero_matrix, mismatch_score):
	# sets up dynamic program for local alignment to fill matrix
	# python implementation of Smith Waterman Algorithm for local alignments

	#print("dp local matrix for {0}".format((sequence_large_genome, sequence_small_align)))
	
	# top is the top of the matrix, and side is the right hand side (location for both sequences)
	top_sequence = list(sequence_large_genome)
	side_sequence = list(sequence_small_align) # breaks sequence into a list of values

	column_value = len(top_sequence)
	#print(column_value == len(zero_matrix[0]))
	row_value = len(side_sequence)
	#print(row_value == len(zero_matrix))

	update_matrix = list(zero_matrix) # re-creates matrix so it can be updated
	max_dict = {} # stores the highest value found in the sequence (stores potential duplicates)

	traceback_dict = {} # {(row, column): ('top/left/diagonal', (row, column previous))}
	location_value_dict = {} # stores the column/row for each value
	
	 # include the top row and side row (i,0), (j,0) -> stopping condition
	for row in range(row_value):
		location_value_dict[(row,0)] = 0
		traceback_dict[(row,0)] = ('intial', (0,0))
	for column in range(column_value):
		location_value_dict[(0,column)] = 0
		traceback_dict[(0,column)] = ('intial', (0,0))

	for row in range(1, row_value): # ignores the first column's dash
		for column in range(1, column_value): # plus one to ignore the first row's dash
			#print("row={0}, column{1}".format(row, column))
			# match or mismatch for current cell
			bases_match = False
			#print("{0}{1}".format(top_sequence[column], side_sequence[row]))
			if top_sequence[column] == side_sequence[row]:
				current_score_index = match_score
				bases_match = True
			else:
				current_score_index = mismatch_score

			# determine the adjacent cells that will be compared
			top = update_matrix[row-1][column]
			left = update_matrix[row][column-1]
			diagonal = update_matrix[row-1][column-1]
			#print("top: {0}+{1}, left: {2}+{1}, diagonal: {3}+{4}".format(top, gap_score, left, diagonal, current_score_index))

			# update adjacent cells if the current index matches
			diagonal += current_score_index
			top += gap_score
			left += gap_score

			#print("[{0},{1}], max = {2}\n".format(row, column, max(0, top, left, diagonal))
			update_matrix[row][column] = max(0, top, left, diagonal)
			
			# updates max_dict for all values, will remove the smallest after filled
			location_value_dict[(row,column)] = update_matrix[row][column]

			# store the path that the cells were populated
			if diagonal >= left and diagonal >= top:
				traceback_dict[(row,column)] = ('diagonal', (row-1, column-1))
			if left >= top and left > diagonal:
				traceback_dict[(row,column)] = ('left', (row, column-1))
			if top > left and top > diagonal:
				traceback_dict[(row,column)] = ('top', (row-1, column))

	# only return the location for the max values (allows for duplicates)
	largest_value_in_matrix = max(location_value_dict.iteritems(), key=operator.itemgetter(1))[1]
	#print(largest_value_in_matrix)
	# can be multiple cells that contain the same value (multiple paths)
	# example: {(3, 9): 9, (3, 12): 9}
	max_dict = { k:v for k, v in location_value_dict.items() if v == largest_value_in_matrix }

	return (update_matrix, traceback_dict, max_dict, location_value_dict)

def traceBackPath(traceback_dictionary, max_location_dictionary, location_value_dictionary):
	# returns a dictionary with the path and starting location for all values
	traceback_word_paths = {} # list of lists for all paths
	internal_path = [] # one path that will be included in the traceback above


	for key in max_location_dictionary:
		#print("\n\t\tKEY: {0}\n".format(key))
		current_location = key
		internal_path.append(traceback_dictionary[current_location][0])
		while (location_value_dict[current_location] != 0):
			#print(current_location)
			#print(traceback_dictionary[current_location])
			#print(location_value_dict[current_location])
			current_location = traceback_dictionary[current_location][1]
			if (traceback_dictionary[current_location][0] != 'intial'): # include full path, ignore the final 0
				internal_path.append(traceback_dictionary[current_location][0])
		traceback_word_paths[key] = internal_path
		internal_path = []

	#print(traceback_word_paths)
	return traceback_word_paths

def alignSequencesStrings(directions, genome_sequence, aligner_sequence):
	# directions are taken in from top left to bottom right (in order of the sequence)
	aligned_small_sequence = ''
	aligned_large_genome = ''

	total_aligned = {} # contains a list of lists for each aligned sequence

	# alignments are populated in reverse (turns given list into a string and reverses)
	# strings are reversed since the directions are in the reverse (populated from the end to the beginning of the string)
	genome_string = (''.join(genome_sequence))[::-1][:-1]
	aligner_string = (''.join(aligner_sequence))[::-1][:-1] #removes the gap unneeded preceding gap

	for key in directions:
		# start location determined by the max value
		# position realigned for the reversed string used by the directions
		position_align = (len(aligner_sequence)-key[0])-1 # minus one accounts for location in row/columns that start at 0
		position_genome = (len(genome_sequence)-key[1])-1
		align_path = directions[key]

		for i in range(len(align_path)):
			# dynamic programming rules for the direction that a path takes
			if align_path[i] == 'diagonal':
				if (position_genome < len(genome_string)) or (position_align < len(aligner_string)):
					aligned_large_genome += genome_string[position_genome]
					aligned_small_sequence += aligner_string[position_align]
					position_genome += 1
					position_align += 1
			if align_path[i] == 'top':
				if (position_align < len(aligner_string)):
					aligned_large_genome += '-'
					aligned_small_sequence += aligner_string[position_align]
					position_align += 1
			if align_path[i] == 'left':
				if (position_genome < len(genome_string)):
					aligned_large_genome += genome_string[position_genome]
					aligned_small_sequence += '-'
					position_genome += 1
		
		total_aligned[key] = (aligned_small_sequence[::-1],aligned_large_genome[::-1]) # reverse strings to return to normal reading
		#print("small: {0}, large: {1}".format(aligned_small_sequence[::-1], aligned_large_genome[::-1]))
		aligned_small_sequence = ''
		aligned_large_genome = ''
	#print(total_aligned)
	return total_aligned
########################################################################
## PRINT ALIGNMENT, MATRIX, etc...
def neatPrintMatrix(matrix, top_sequence, side_sequence):
	# print with character inline, top is genome, side is sequence to be aligned
	'''
	dp local matrix for ('-TATAGACACATACG', '-CAT')
		   -  T  A  T  A  G  A  C  A  C  A  T  A  C  G
		- [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		C [0, 0, 0, 0, 0, 0, 0, 3, 2, 3, 2, 1, 0, 3, 2]
		A [0, 0, 3, 2, 3, 2, 3, 2, 6, 5, 6, 5, 4, 3, 2]
		T [0, 3, 2, 6, 5, 4, 3, 2, 5, 4, 5, 9, 8, 7, 6]
	'''
	print("\n   {0}".format("  ".join(top_sequence)))
	for i in range(len(matrix)):
		print("{0} {1}".format(side_sequence[i], matrix[i]))

def findRangeConfidence(high_align_dictionary, low_align_dictionary, genome_sequence, aligner_sequence):
	# find range of values that contain the sequence in a variety of confidences
	# high confidence: exact match is found (could be found in multiple locations across the genome for a larger range than sequence length)
	# low confidence: no exact matches were found (range includes mismatches)
	dict_min_range_low = min([int(i[1]) for i in high_align_dictionary.keys()])
	dict_max_range_low = max([int(i[1]) for i in high_align_dictionary.keys()])
	dict_seq_length = max([len(i[1]) for i in high_align_dictionary.values()])

	if high_align_dictionary == low_align_dictionary:
		match_es = 'match'
		if (len(high_align_dictionary.values()) != 1):
			match_es = 'matches'
		print("exact {0} found").format(match_es)
		print("HIGH CONFIDENCE RANGE: Between {0} to {1}".format(dict_min_range_low-dict_seq_length+1, dict_max_range_low))
	else:
		min_range_low = min([int(i[1]) for i in low_align_dictionary.keys()])
		min_range_high = min([int(i[1]) for i in high_align_dictionary.keys()])
		min_range = min(min_range_low, min_range_high)
		
		max_range_low = max([int(i[1]) for i in low_align_dictionary.keys()])
		max_range_high = max([int(i[1]) for i in high_align_dictionary.keys()])
		max_range = max(max_range_low, max_range_high)
		
		seq_length_low = max([len(i[1]) for i in low_align_dictionary.values()])
		seq_length_high = max([len(i[1]) for i in high_align_dictionary.values()])
		max_seq_length = max(seq_length_low, seq_length_high)
		print("no exact matches found")
		print("HIGH CONFIDENCE RANGE: Between {0} to {1}".format(dict_min_range_low-dict_seq_length+1, dict_max_range_low))
		print("LOW  CONFIDENCE RANGE: Between {0} to {1}".format(max_seq_length-min_range, max_range))
	print("\n")

def printNeatAlignment(high_align_dictionary, low_align_dictionary, genome_sequence, aligner_sequence, genome_filename, genome_name, align_name):
	# print alignments for console
	'''
	####################################################################
	Genome File: <genome_filename>
	Genome Seq:  <genome_name>
	Align Seq:   <aligner_sequence_name>

	Genome Seq Length: <length of genome>
	Total matches found: <total times the sequence appears>
	
	<genome_name> csci4314 match <percent match> ID=<#> + <start> <end>
			CAT-AT--GGGA
			|||~||~~|
	  TATAGACATCATACG
	...
	####################################################################
	EXACT MATCH FOUND EXAMPLE:
	Genome File: genome_example.fasta
	Genome Seq:  GEN_0
	Align Seq:   SEQ_2

	Genome Seq Length: 14
	Total matches found: 1

	exact match found
	HIGH CONFIDENCE RANGE: Between 9 to 11


	GEN_0  csci4314  match  100.00  ID=1  +  9  11

	Full Sequence Display:
			  CAT
			  |||
	  TATAGACACATACG
	####################################################################
	NO EXACT MATCH FOUND EXAMPLE:
	Genome File: genome_example.fasta
	Genome Seq:  GEN_0
	Align Seq:   SEQ_1

	Genome Seq Length: 14
	Total matches found: 3

	no exact match found
	HIGH CONFIDENCE RANGE: Between 0 to 10
	LOW  CONFIDENCE RANGE: Between 2 to 13


	GEN_0  csci4314  match  100.00  ID=1  +  2  10

	Full Sequence Display:
		  CGCG-TA-A-A-AAAAA
		  ~||~|~|~|
	  TATAGACACATACG

	GEN_0  csci4314  match  30.00  ID=2  +  4  13

	Full Sequence Display:
		 CGCGCGCGTAAAA
		 XXX|XX|X|X
	  TATAGACACATACG

	GEN_0  csci4314  match  100.00  ID=3  +  1  8

	Full Sequence Display:
		   CGCGCTA-A-A-AAAA
		   ||~|~|~|
	  TATAGACACATACG
	'''
	# combine low/high confidence dictionaries into one that can be iterated through
	total_updated_dictionary = dict(high_align_dictionary)
	total_updated_dictionary.update(low_align_dictionary)
	print("####################################################################")
	print("Genome File: {0}".format(genome_filename))
	print("Genome Seq:  {0}".format(genome_name))
	print("Align Seq:   {0}".format(align_name))
	print("\nGenome Seq Length: {0}".format(len(genome_sequence)-1))
	print("Total matches found: {0}\n".format(len(total_updated_dictionary.keys()))) # total number of times a sequence appears
	
	findRangeConfidence(high_align_dictionary, low_align_dictionary, genome_sequence, aligner_sequence) # print range of confidences
	
	aligner_sequence = aligner_sequence[1:] # remove preceding gapping done for sequence alignment
	genome_sequence = genome_sequence[1:]

	id_counter = 1 # keeps track for printing the id value for each found value

	for key in total_updated_dictionary:
		both_alignments = total_updated_dictionary[key]
		small_align = both_alignments[0]
		large_genome = both_alignments[1]

		print(large_genome)
		symbols = symbolGenerate(large_genome, small_align)
	
		percent_match = float(len(large_genome)-symbols.count('X'))/float(len(large_genome))*100 # counts mismatches
		#print("Starting Location: {0}".format(key))
		#print("{0} - {1} = {2}/{3} = {4}".format(len(large_genome), symbols.count('X'), len(large_genome)-symbols.count('X'), len(large_genome), percent_match))
		print("{0}  csci4314  match  {1:.2f}  ID={2}  +  {3}  {4}".format(genome_name, percent_match, id_counter,  key[1]-len(large_genome)+1, key[1]))
		id_counter += 1
		
		
		# commented out for simple printing
		# print with only aligned section
		#printSequenceAligned(small_align, large_genome, symbols)
		# print with sequence displayed
		printFullSequence(key, small_align, aligner_sequence, large_genome, genome_sequence, symbols)
		
	print("####################################################################")

def printBowtieFormat(high_align_dictionary, align_name):
	# print format to match the output of bowtie to be spliced
	# if there is more than one match for the highest value
	if len(high_align_dictionary) > 1:
		multiple_start_locations = high_align_dictionary.keys()
		min_loc = min(multiple_start_locations, key=operator.itemgetter(1))
		small_align = high_align_dictionary[min_loc][0]
		large_genome = high_align_dictionary[min_loc][1]
		
		genome_start_position = max(0, min_loc[1] - len(large_genome))+1
		genome_end_position = min_loc[1]
	# if there is only one match for the highest value
	else:
		for key in high_align_dictionary:
			small_align = high_align_dictionary[key][0]
			large_genome = high_align_dictionary[key][1]
			
			genome_start_position = max(0, key[1] - len(large_genome))+1
			genome_end_position = key[1]
	
	total_gaps = large_genome.count('-')
	symbols =  symbolGenerate(large_genome, small_align)
	percent_match = float(len(large_genome)-symbols.count('X'))/float(len(large_genome))*100 # counts mismatches
	print("{0}\t+\tSTART={1}\tEND={2}\tGAPS={3}\tMATCH={4}".format(align_name, genome_start_position, genome_end_position, total_gaps, percent_match))

def symbolGenerate(large_genome, small_align):
		symbols = ""
		# set up symbols between both alignments when printed
		for i in range(len(large_genome)):
			if (small_align[i] == large_genome[i]):
				symbols = symbols + '|'
			elif (small_align[i] == '-' or large_genome[i] == '-'):
				symbols = symbols + '~'
			elif (small_align[i] != '-' and large_genome[i] != '-') and (small_align[i] != large_genome[i]):
				symbols = symbols + 'X'
		return symbols

def printSequenceAligned(small_align, large_genome, symbols):
		print("\nAlignment Sequence Only:")
		print("small_align_: {0}".format(small_align))
		print("              {0}".format(symbols))
		print("large_genome: {0}".format(large_genome))
		'''
		Example:
		small_align_: CAT-AT--G
              |||~||~~|
		large_genome: CATCATACG
		When the full looks like: 'TATAGACATCATACG', 'CATATGGGA'
		'''

def printFullSequence(key, small_align, aligner_sequence, large_genome, genome_sequence, symbols):
		#### PRINT FULL SEQUENCE
		aligned_start_position = max(0, key[0] - len(small_align)) # produces a non-negative range
		aligned_end_position = aligned_start_position + len(small_align.replace('-', ''))
		#print(range(aligned_start_position, aligned_end_position))

		updated_aligner_sequence = ''
		updated_aligner_sequence = aligner_sequence[0:aligned_start_position] + small_align + aligner_sequence[aligned_end_position:]
		#updated_aligner_sequence = aligner_sequence[0:aligned_start_position] + '^' + small_align + '^' + aligner_sequence[aligned_end_position:]
		# uses '^' to seperate sequences for easier reading
		#print("updated_aligner: {0}\n".format(updated_aligner_sequence))
		
		genome_start_position = max(0, key[1] - len(large_genome))
		genome_end_position = genome_start_position + len(large_genome.replace('-', ''))
		#print(range(genome_start_position, genome_end_position))

		updated_genome_sequence = ''
		updated_genome_sequence = genome_sequence[0:genome_start_position] + large_genome + genome_sequence[genome_end_position:]
		#updated_genome_sequence = genome_sequence[0:genome_start_position] + '^' + large_genome + '^' + genome_sequence[genome_end_position:]
		#print("updated_genome_: {0}".format(updated_genome_sequence))

		spaces =  ' '*max(aligned_start_position, genome_start_position) # + ' ' # aligns the empty spaces to line up sequences on print console
		# plus one accounts for the use of '^' in the print statement if used


		print("\nFull Sequence Display:")
		print("  {0}{1}".format(spaces, updated_aligner_sequence))
		print("  {0}{1}".format(spaces, symbols))
		print("  {0}\n".format(updated_genome_sequence))
		'''
		Example:
		Full Sequence Display:
		CATA-TGGGA
		X|||~XXXX|
		TATAGACACATACG
		'''

########################################################################

if __name__ == '__main__':
	import argparse # parse given flag options
	parser = argparse.ArgumentParser(description="flag format given as: -SF <sequence_filename> -S <sequence> -AF <align_filename> -A <align sequence>")
	parser.add_argument('-SF', '-Sequence_File', help="given sequence file, must be .fasta")
	parser.add_argument('-S', '-Sequence_Name', help="given sequence file found in sequence file")
	parser.add_argument('-AF', '-Alignment_File', help="given sequence to align's file, must be .fasta")
	parser.add_argument('-P', '-Print_Format', help="how to print the final results")
	# order of arguments does not matter

	args = parser.parse_args()
	sequence_filename = args.SF
	sequence_name = args.S
	alignment_filename = args.AF
	print_format = args.P
	
	arguments = [sequence_filename, sequence_name, alignment_filename, print_format]
	# if either arguments are not given (left empty), then quit
	if None in arguments:
		if sequence_filename is None:
			print("sequence source filename not given")
			exit()
		if sequence_name is None:
			print("sequence_name not given")
			exit()
		if alignment_filename is None:
			print("alignment source filename not given")
			exit()
		if print_format is None:
			print("print format is not provided")
			exit()

	# only accept file that are fasta
	if sequence_filename.endswith('.fasta') and alignment_filename.endswith('fasta'):
		pass
	elif sequence_filename.endswith('.FASTA') and alignment_filename.endswith('.FASTA'):
		pass
	elif sequence_filename.endswith('.FA') and alignment_filename.endswith('.FA'):
		pass
	elif sequence_filename.endswith('.fa') and alignment_filename.endswith('.fa'):
		pass
	else:
		if not sequence_filename.endswith('.fasta') or sequence_filename.endswith('.FASTA') or sequence_filename.endswith('.fa'):
			print("\n\t{0} is not a .fasta file, please choose a different file\n".format(sequence_filename))
			exit()
		if not alignment_filename.endswith('.fasta') or alignment_filename.endswith('.FASTA') or alignment_filename.endswith('.fa'):
			print("\n\t{0} is not a .fasta file, please choose a different file\n".format(alignment_filename))
			exit()
	'''
	The following returns a list of the sequence headers ['chrI', 'chrII',
	etc...) and genome ['ATC', 'TGGC', etc..] that is spliced out of a
	given file based on if the line starts with > as seen in fasta
	'''

	small_align_dict = readingFileDict(alignment_filename)
	#print(small_align_dict)
	large_genome_dict = readingFileDict(sequence_filename)
	#print(large_genome_dict)
	# returns a dictionary that contains the name and sequence that will be compared against with all the genomes to be compared against
	
	# exit if the genome name is not found in the genome file
	if sequence_name not in large_genome_dict.keys():
		print("\n\t'{0}' is not found in genome file '{1}', choose a different file or sequence name".format(sequence_name, sequence_filename))
		print("\n\tAvailable genomes names:")
		for key in large_genome_dict.keys():
			print("\t{0}".format(key))
		print("\n")
		exit()

	paired_seq = orderedPairs(small_align_dict, sequence_name) # creates pairs for each sequence and the genome used

	# creates a zero matrix for each aligned sequences compared to the larger genome
	for pair in paired_seq:
		#print(pair)
		aligned_sequence = small_align_dict[pair[0]]
		genome_to_align = large_genome_dict[pair[1]]
		#print("aligned_sequence: {0}, genome: {1}".format(aligned_sequence, genome_to_align))
		
		# if the smaller aligning sequence is greater than the genome, exit
		if len(aligned_sequence) > len(genome_to_align):
				print("\n\t'{0}' is larger than the genome '{1}' to compare against, choose a different sequence or genome\n".format(pair[0], pair[1]))
				exit()

		# add gap value at the start of both sequences
		aligned_sequence = '-' + aligned_sequence
		genome_to_align = '-' + genome_to_align
		# width of matrix is the size of the genome, height is the size of the aligner sequence
		zero_matrix = zeroMatrix(len(genome_to_align), len(aligned_sequence))

		high_confidence_mismatch_score = -2 # limits the amount of mismatches

		dp_total = DPlocalMatrix(genome_to_align, aligned_sequence, zero_matrix, high_confidence_mismatch_score)
		
		dp_local_matrix = dp_total[0] # the matrix produced by the dynamic program
		#neatPrintMatrix(dp_local_matrix, genome_to_align, aligned_sequence)

		traceback_dict = dp_total[1] # contains the path that populated the values
		max_location_dict = dp_total[2] # contains a dictionary that stores the row/column of the largest value in the matrix
		location_value_dict = dp_total[3] # contains a dictionary that stores the row/column and it's associated value
		traceback_path = traceBackPath(traceback_dict, max_location_dict, location_value_dict)
		high_confidence_alignment_dicts = alignSequencesStrings(traceback_path, genome_to_align, aligned_sequence)

		# set up the range for high and low confidence based on found values
		low_confidence_mismatch_score = 2 # expand total mismatches
		low_dp_total = DPlocalMatrix(genome_to_align, aligned_sequence, zero_matrix, low_confidence_mismatch_score)

		low_dp_local_matrix = low_dp_total[0] # the matrix produced by the dynamic program
		low_traceback_dict = low_dp_total[1] # contains the path that populated the values
		low_max_location_dict = low_dp_total[2] # contains a dictionary that stores the row/column of the largest value in the matrix
		low_location_value_dict = low_dp_total[3] # contains a dictionary that stores the row/column and it's associated value

		low_traceback_path = traceBackPath(low_traceback_dict, low_max_location_dict, low_location_value_dict)
		low_confidence_alignment_dicts = alignSequencesStrings(low_traceback_path, genome_to_align, aligned_sequence)
		
		if print_format.upper() == 'USER':
			# user friendly output:
			printNeatAlignment(high_confidence_alignment_dicts, low_confidence_alignment_dicts, genome_to_align, aligned_sequence, sequence_filename, pair[1], pair[0])
		if print_format.upper() == 'BOWTIE':
			# bowtie-format output:
			printBowtieFormat(high_confidence_alignment_dicts, pair[0])
