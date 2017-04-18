import os
import re
directory = "../test_results"

# number of reads per file
TR = 50
#specify standard deviation for accepted alignments
SD = .1

# loop through each file in directory
for filename in os.listdir(directory):
    #counter for successfully matches
    hits = 0
    false_positives = 0
    true_positives = 0
    false_negatives = 0
    if filename.endswith(".fq"):
        rd = int(re.search('_rd(.+?).fq', filename).group(1))
        working_file = open(filename, 'r')
        # loop through each file
        for line in working_file:

            #bowtie found a match
            hits+=1
            elements = line.split("\t")
            bowtie_strand = elements[1]
            bowtie_position = int(elements[3])

            title_info = elements[0].split()
            real_alignment = title_info[2]
            strand = ""
            alignment = 0

            # regex search for actual alignment if complement
            x = re.search('position=complement\((.+?)\.\.', real_alignment)
            if x:
                strand = '-'
                alignment = int(x.group(1))
            else:
                 strand = '+'
                 # strand is not complement so regex search for + strand
                 alignment = int(re.search('position=(.+?)\.\.',real_alignment).group(1))

            #calculate allowed mis-alignment
            shift_down = alignment - (rd*SD)
            shift_up = alignment + (rd*SD)
            if bowtie_strand == strand and bowtie_position >= shift_down and bowtie_position <= shift_up:
                # the alignment bowtie choose was right!
                true_positives += 1
            else:
                # the alignment bowtie choose was wrong
                false_positives += 1
                # print "true: " + str(bowtie_position) + " low: " +str(shift_down) + " high: " + str(shift_up) + " predicted: " + str(bowtie_position)
        false_negatives = TR - hits
        print ">" + filename
        print "true positives: " + str(true_positives) + " false positives: " + str(false_positives) + " false negatives: " + str(false_negatives)
