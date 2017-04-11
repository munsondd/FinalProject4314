########################################################################
# Final Project: CSCI-4314
    # Date: 4/10/2017
	# Benchmark sequence aligner against Bowtie

    # Format:
        # File contain sequence, sequence name
        	# -F <fasta_filename> -S <sequence_name>
########################################################################

if __name__ == '__main__':
	import argparse # parse given flag options
	parser = argparse.ArgumentParser(description="flag format given as: -F <filename> -L <length/size of kmer")
	parser.add_argument('-F', '-filename', help="given filename with sequence, must be .fasta")
	parser.add_argument('-S', '-Sequence', help="given sequence name")
	# order of arguments does not matter

	args = parser.parse_args()
	filename = args.F
	sequence_name = args.S

	arguments = [filename, sequence_name]
	# if either arguments are not given (left empty), then quit
	if None in arguments:
		if filename is None:
			print("filename not given")
			exit()
		if sequence_name is None:
			print("sequence_name not given")
			exit()

	# only accept file that are fasta
	if filename.endswith('.fasta'):
		pass
	elif filename.endswith('.FASTA'):
		pass
	else:
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
	with open(filename, "r") as given_file:
		seq = '' # helped by Subi Nair
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

	seq_dict = seqGenDict(header_list, nucleo_list) # returns a dictionary that combines header/genome {seq:gen}
	# O(N) and size 2n

	if all(i is '[^ATCG]' for i in nucleo_list): # exit if nucleotides contain non-ATCG characters (not A, T, C or G)
		print("\n\tSequence contains non-nucleotide character, therefore, no matches found\n")
		exit()
