########################################################################
# Final Project: CSCI-4314
    # Date: April 2017
    
    # Purpose: Benchmark sequence aligner against Bowtie
    # Create a local aligner to benchmark against Bowtie
    # Compares a list of small sequences to one large genome and returns all alignments
    # Finds local regions with the highest levels of conservation/similarity

    # Psuedo-Code:
    
    # Format: -SF <sequence_filename> -S <sequence> -AF <align_filename>
	# Sequence file: File that contains all large genomes that will be aligned against
		# Sequence: The specific sequence within the file that will be aligned against
	# Alignment File: File that contains all the smaller sequences to be aligned
########################################################################

import argparse # parse given flag options
import re
'''
scoring_matrix = {'AA':  0, 'AC': 1, 'AG':  1, 'AT': 1, 'A-': 1,
				 'CA': 1, 'CC':  0, 'CG': 1, 'CT':  1, 'C-': 1,
				 'GA':  1, 'GC': 1, 'GG':  0, 'GT': 1, 'G-': 1,
				 'TA': 1, 'TC':  1, 'TG': 1, 'TT':  0, 'T-': 1,
				 '-A': 1, '-C': 1, '-G': 1, '-T': 1, '--': 1}
'''
match_score = 3
mismatch_score = -2
gap_score = -1

'''
scoring_matrix = {'AA':  0, 'AC': 1, 'AG':  1, 'AT': 1, 'A-': 1,
				 'CA': 1, 'CC':  0, 'CG': 1, 'CT':  1, 'C-': 1,
				 'GA':  1, 'GC': 1, 'GG':  0, 'GT': 1, 'G-': 1,
				 'TA': 1, 'TC':  1, 'TG': 1, 'TT':  0, 'T-': 1,
				 '-A': 1, '-C': 1, '-G': 1, '-T': 1, '--': 1}
'''

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

def DPlocalMatrix(sequence_large_genome, sequence_small_align, zero_matrix):
	# sets up dynamic program for local alignment to fill matrix
	print("dp local matrix for r {0}".format((sequence_large_genome, sequence_small_align)))
	
	# top is the top of the matrix, and side is the right hand side (location for both sequences)
	top_sequence = list(sequence_large_genome)
	side_sequence = list(sequence_small_align) # breaks sequence into a list of values

	column_value = len(top_sequence)
	#print(column_value == len(zero_matrix[0]))
	row_value = len(side_sequence)
	#print(row_value == len(zero_matrix))

	update_matrix = list(zero_matrix) # re-creates matrix so it can be updated

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
	return update_matrix

def smithWaterman(sequence_large_genome, sequence_small_align, dp_matrix):
	# python implementation of Smith Waterman Algorithm for local alignments
	print("smith waterman for {0}".format((sequence_large_genome, sequence_small_align)))

########################################################################
## PRINT ALIGNMENT, MATRIX, etc...
def neatPrint(matrix, top_sequence, side_sequence):
	print("\n   {0}".format("  ".join(top_sequence)))
	for i in range(len(matrix)):
		print("{0} {1}".format(side_sequence[i], matrix[i]))

########################################################################

if __name__ == '__main__':
	import argparse # parse given flag options
	parser = argparse.ArgumentParser(description="flag format given as: -SF <sequence_filename> -S <sequence> -AF <align_filename> -A <align sequence>")
	parser.add_argument('-SF', '-Sequence_File', help="given sequence file, must be .fasta")
	parser.add_argument('-S', '-Sequence_Name', help="given sequence file found in sequence file")
	parser.add_argument('-AF', '-Alignment_File', help="given sequence to align's file, must be .fasta")
	# order of arguments does not matter

	args = parser.parse_args()
	sequence_filename = args.SF
	sequence_name = args.S
	alignment_filename = args.AF
	
	arguments = [sequence_filename, sequence_name, alignment_filename]
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

	# only accept file that are fasta
	if sequence_filename.endswith('.fasta') and alignment_filename.endswith('fasta'):
		pass
	elif sequence_filename.endswith('.FASTA') and alignment_filename.endswith('.FASTA'):
		pass
	else:
		if not sequence_filename.endswith('.fasta') or sequence_filename.endswith('.FASTA'):
			print("\n\t{0} is not a .fasta file, please choose a different file\n".format(filename))
			exit()
		if not alignment_filename.endswith('.fasta') or alignment_filename.endswith('.FASTA'):
			print("\n\t{0} is not a .fasta file, please choose a different file\n".format(filename))
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
		aligned_sequence = small_align_dict[pair[0]]
		genome_to_align = large_genome_dict[pair[1]]
		
		# if the smaller aligning sequence is greater than the genome, exit
		if len(aligned_sequence) > len(genome_to_align):
				print("\n\t'{0}' is larger than the genome '{1}' to compare against, choose a different sequence or genome\n".format(pair[0], pair[1]))
				exit()

		# add gap value at the start of both sequences
		aligned_sequence = '-' + aligned_sequence
		genome_to_align = '-' + genome_to_align

		# width of matrix is the size of the genome, height is the size of the aligner sequence
		zero_matrix = zeroMatrix(len(genome_to_align), len(aligned_sequence))
		print(pair)
		#for row in zero_matrix:
		#	print(row)
		dp_local_matrix = DPlocalMatrix(genome_to_align, aligned_sequence, zero_matrix)
		neatPrint(dp_local_matrix, genome_to_align, aligned_sequence)
		#smithWaterman(genome_to_align, aligned_sequence, dp_local_matrix)
		print("\n")
