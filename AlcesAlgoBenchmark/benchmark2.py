import os
import re
directory = "../AlcesAlgoBenchmark"

# number of reads per file
TR = 50
#specify standard deviation for accepted alignments
SD = .1
tableset = []
# loop through each file in directory
for filename in os.listdir(directory):
    #counter for successfully matches
    hits = 0
    false_positives = 0
    true_positives = 0
    false_negatives = 0
    if filename.endswith(".txt"):
        # rd is now specified as len?
        rd = int(re.search('_len_(.+?)_r', filename).group(1))
        working_file = open(filename, 'r')
        # loop through each file
        for line in working_file:

            #found a match
            hits+=1
            elements = line.split("\t")
            our_strand = elements[1]

            start_position = int(re.search('START=(.*)',elements[2]).group(1))
            end_position = int(re.search('END=(.*)',elements[3]).group(1))

            title = elements[0]
            strand = ""
            left_alignment = 0
            right_alignment = 0

            # regex search for actual alignment if complement
            x = re.search('POSITION=COMPLEMENT\((.+?)\)DESCRIPTION', title)
            if x:
                strand = '-'
                alignment = x.group(1)
                left_alignment = int(re.search('(.+?)\.\.', alignment).group(1))
                right_alignment = int(re.search('\.\.(.+?)', alignment).group(1))
            else:
                 strand = '+'
                 # strand is not complement so regex search for + strand
                 alignment = re.search('POSITION=(.+?)DESCRIPTION', title).group(1)
                 left_alignment = int(re.search('(.+?)\.\.', alignment).group(1))
                 right_alignment = int(re.search('\.\.(.*)', alignment).group(1))

            #calculate allowed mis-alignment
            #print "our strand: " + str(our_strand) + " real strand: " + str(strand)
            #print "our start: " + str(start_position) + " real start: " + str(left_alignment)
            #print "our end: " + str(end_position) + " real end: " + str(right_alignment)
            if our_strand == strand and left_alignment == start_position:
                # the alignment was right!
                true_positives += 1
            else:
                # the alignment was wrong
                false_positives += 1
                # print "true: " + str(bowtie_position) + " low: " +str(shift_down) + " high: " + str(shift_up) + " predicted: " + str(bowtie_position)
        false_negatives = TR - hits

        #calculating precision and recall
        precision = 0
        recall = 0
        if true_positives + false_positives != 0:
            precision = float(true_positives)/float((true_positives + false_positives))
        if true_positives + false_negatives != 0:
            recall = float(true_positives)/float((true_positives + false_negatives))
        #calculating f_value
        f_value = 0
        if recall + precision != 0:
            f_value = (2* recall * precision)/ (recall + precision)
        tableset.append((rd,true_positives*2,false_positives*2,false_negatives*2,f_value))

        print ">" + filename
        print "true positives: " + str(true_positives) + " false positives: " + str(false_positives) + " false negatives: " + str(false_negatives)
        print "precision: " + str(precision) + " recall: " + str(recall) + " f value: " + str(f_value)
tableset = sorted(tableset)
print "bp\ttp\tfp\tfn\tf value"
for row in tableset:
    print str(row[0]) + "\t" + str(row[1]) + "\t" + str(row[2]) +"\t" + str(row[3]) + "\t" + str(row[4])
