########################################################################
# Final Project: CSCI-4314
    # Date: April 2017
    
    # Purpose: Benchmark sequence aligner against Bowtie
    # Create a local aligner to benchmark against Bowtie
    # Compares a list of small sequences to one large genome and returns all alignments

    # Psuedo-Code:
    
    # Format: -SF <sequence_filename> -S <sequence> -AF <align_filename> -A <align sequence>
	# Sequence file: File that contains all large genomes that will be aligned against
		# Sequence: The specific sequence within the file that will be aligned against
	# Alignment File: File that contains all the smaller sequences to be aligned
########################################################################

import argparse # parse given flag options
import re

def seqDictPairs(header_list, sequence_list):
	# creates a dictionary between the sequence (header) and the associated genome {seq:genome} dictionary
	# takes in the header list and the sequence list to combine into a single diectionary that can be searched
	# setting up sequence genome dictionary
	seq_gen_dict = {}
	seq_gen_dict = zip(header_list, sequence_list) # combine the two lists
	seq_gen_dict = dict(seq_gen_dict) # create new dictionary from the lists
	return seq_gen_dict

def SmithWaterman(sequence_large_genome, sequence_small_align):
	# python implementation of Smith Waterman Algorithm for local alignments
	pass


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

	fullList = []
	seqlist = []
	header_list = []
	nucleo_list = []

	append = fullList.append # avoid re-using append (to improve running time)
	with open(alignment_filename, "r") as given_file:
		seq = ''
		for line in given_file:
			line = line.rstrip('\n').replace(" ", "@").replace("\t", "@") # replace spaces with known character and replace tabs
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
				element = element.replace("@", " ")
				element = element.strip(">") # assumes all header/sequences starts with >
				append(element) # returns a list of  headers ['chrI', 'chrII', etc...]

		# Returns the genome/sequence
		nucleo_list =  [x for x in seqList if '>' not in x] # returns a list of sequence ['ATC', 'TGGC', etc..]
		nucleo_list = map(str.upper, nucleo_list) # convert all sequences to upper case for consitency

	small_align_dict = seqDictPairs(header_list, nucleo_list)
	# returns a dictionary that combines header/align {seq:alg} for the file with all the short sequences that will be aligned
	#print(small_align_dict)

	# find the sequence that is requested in the arguments, if not found, exit
	found = False
	with open(sequence_filename, "r") as given_file:
		sequence_name = sequence_name.lower()
		if not found:
			for line in given_file:
				if line.startswith('>'):
					search_line = re.search((r'{0}'.format(sequence_name)), line.lower())
					if search_line is not None:
						found=True
						sequence_in_file = next(given_file) # returns the next line in the file which is the associated sequence
		if not found: # if the value is never found, exit
			print("\n\t{0} is not found in given sequence file {1}, please choose a different file or sequence name\n".format(sequence_name, sequence_filename))
			exit()
	large_genome_dict = {}
	sequence_in_file = sequence_in_file.strip('\n')
	large_genome_dict[sequence_name] = sequence_in_file
	# returns a dictionary that contains the name and sequence that will be compared against with all the genomes to be compared against
	#print(large_genome_dict)
