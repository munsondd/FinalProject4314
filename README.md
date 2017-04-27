# Breakdown of Arguments

```
Format: -SF <genome filename> -S <specific genome> -AF <align filename> -P <print format>
```
#### -SF: Genome Filename
The genome file is a .fasta file that contains a list of different genomes. If the file only contains one genome, it will still need to include both a filename and a specific genome as listed below in -S
#### -S: Specific Genome
The specific genome within the larger -SF (genome file) that the user wants to compare all the sequences against
#### -AF: Sequence Aligners Filename
A .fasta file that contains all the smaller sequences that the user wants to align against the larger genome
#### -P: Print Format
The specific print formatting that the user wants to output to the console. Currently works for 'user' and 'bowtie'. User is a user friendly printing style and bowtie is the bowtie format used for benchmarking. See below for examples of both


# Run Existing Algorithm
 To run code:
 1. Download repo
 2. Requires python:  
 ```python csci_4314_algo.py -SF genome_example.fasta -S GEN_0 -AF aligners_example.fasta```
 #### Example Output if found EXACTLY (single/multiple instances)
 ```
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

Genome File: genome_example.fasta
Genome Seq:  GEN_0
Align Seq:   SEQ_4

Genome Seq Length: 14
Total matches found: 3

exact matches found
HIGH CONFIDENCE RANGE: Between 1 to 12


GEN_0  csci4314  match  100.00  ID=1  +  11  12

Full Sequence Display:
            TA
            ||
  TATAGACACATACG

GEN_0  csci4314  match  100.00  ID=2  +  1  2

Full Sequence Display:
  TA
  ||
  TATAGACACATACG

GEN_0  csci4314  match  100.00  ID=3  +  3  4

Full Sequence Display:
    TA
    ||
  TATAGACACATACG

 ```
  #### Example Output if found NOT EXACTLY
```
Genome File: genome_example.fasta
Genome Seq:  GEN_0
Align Seq:   SEQ_0

Genome Seq Length: 14
Total matches found: 3

no exact matches found
HIGH CONFIDENCE RANGE: Between 0 to 8
LOW  CONFIDENCE RANGE: Between 3 to 12


GEN_0  csci4314  match  44.44  ID=1  +  2  10

Full Sequence Display:
   ATGAGGAGA
   ||XXXX|X|
  TATAGACACATACG

GEN_0  csci4314  match  88.89  ID=2  +  0  8

Full Sequence Display:
  ATGAGGAGA
  ||~|~||X|
  AT-A-GACAACATACG

GEN_0  csci4314  match  72.73  ID=3  +  2  12

Full Sequence Display:
   AT-GA-GGAGA
   ||~||~XX|X|
  TATAGACACATACG
```
 
 Legend: 
 ```
 '|' = match
 '~' = gap
 'X' = mismatch
 ```
 

### interpret bowtie results

bowtie results are printed as follows:

bowtie results are split into two parts. First is the information about the test file pulled from the fasta title. it is formated as follows:
```<id> reference=<reference info> |position=<expected position> description=<description>```
Next is the bowtie results in the following format:
```<strand{+|-}    <reference info>        <bowtie alignment position>     <sequence>```

# Project Proposal: Benchmark Bowtie

![](http://bowtie-bio.sourceforge.net/images/bowtie_logo.png)


## Summary of Topic
Benchmark sequence alignment algorithm against Bowtie for accuracy with a range of errors in the tested sequence

## The Problem
DNA sequences are made up of hundreds of thousands, millions and even billions of repeated basepairs. The human genome is made up of roughly 3.2 bases of the same repeated four base pairs, ‘A’ ‘T’ ‘G’ ‘C’. Sequence aligners like Bowtie are responsible for finding and aligning a sequence with the reference genome. This can be increasingly difficult for long genomes as well as mutations in the genome. Minor mutations in the genome can lead to false negatives, however to account for all possible values in an alignment is an exponentially difficult task.

## Bowtie
Bowtie is alignment algorithm that takes advantage of the compression ability of the Burrows-Wheeler transformation and FM-index. Bowtie is a sequence aligner that aligns shorter sequences of DNA compared to a large reference genome. Any sequence aligner will encountering an initial bottleneck when attempting to reference or even store a genome of any significant size. The Burrows-Wheeler transform is used in Bowtie for compression and indexing of DNA sequences.

#### Burrows-Wheeler Transformation
Burrows-Wheeler transformation is an algorithm that compresses text. A string of text translated into a compressed reverse permutation of the original. When two of the same letter are aligned, it can be further compressed by tracking the number of times a certain character appears. DNA is made up of only four base pairs, ATGC. When a DNA sequence is compressed via Burrows-Wheeler transformation, the chance that any base pair will be placed next to a the same character character is 1 in 4, or 25% of the time, significantly higher than 1 in 26 or  ~3.84% if Burrows-Wheeler was compressing a sequence of the the same length written in English. As a result of basic probability, Burrows-Wheeler will dramatically reduce the size of a compressed DNA sequence. Sequence aligners like Bowtie and SOAP2 both take advantage of Burrows-Wheeler to reduce memory and storage requirements while still maintaining the original sequence in a reversible format.


#### FM-Index
In order to reverse a Burrows-Wheeler transformation it is important to keep track of an individual character’s rank known as LF Mapping. The rank of a character is equal to the number of times that the character has already appeared. In the sequence “BOB”, the rank of each character would be visually represented by a subscript as “B0O0B1”. The rank doesn’t affect how the sequence is compressed, but it allows for an accurate reversal. The FM-index was published in 2000, to describe how to combine the Burrows-Wheeler transformation and rank to form a index with an efficient use of limited memory. FM-index combines the Burrows-Wheeler transformation with a few additional auxiliary data structures to form an index (Burrows and Wheeler 2000). The FM-index is queried in reverse order based on rank until either the sequence is found or did not exist in the original sequence. If the sequence is found, then the range of the sequence is based on the location of the final characters in the first column of the Burrows-Wheeler transformation.

The dramatic increase in speed that Bowtie demonstrates is the result of combining the aligner with the Burrows-Wheeler transformation that allows for more condensed referencing for improved speed and memory use, beating out other tools which can take hours to perform the same alignment. 

This algorithm first indexes the larger reference genome via Burrows-Wheeler transformation. This dramatically shrinks the memory usage and speed needed to access. For the human genome, the memory footprint is about 1.3 GB at alignment time (Langmead et al. 2009). However, it still needs to account for DNA mismatches like insertions, deletions or mismatches. Currently, to search for a sequence in Burrows-Wheeler it must be an exact match. This algorithm accounts for mismatches via a backtracking search seen in Figure 1. To allow for mismatches “Bowtie conducts a quality-aware, greedy, randomized, depth-first search through a space of possible alignments” ((Langmead et al. 2009). 

To avoid time consuming or any excessive amount of backtracking, the algorithm favors a ‘double indexing’ strategy (Langmead et al. 2009). This is done by creating two indices for the genome, containing the ‘forward’ index and the ‘mirror’ or ‘reverse’ index.
This method has yielded significant results in running time in comparison to two similar aligners (SOAP and Maq). For 36, 50 and 76 bp, Bowtie showed at least a 14.9x increase in running time. For 36 bp, Bowtie improve on the approximately four hour running time of Maq and the 17 hour running time of SOAP with a CPU running time of 6 minutes and 15 seconds (Langmead et al. 2009). Bowtie is not a perfect replacement for other aligners like BLAST however. Bowtie is designed for short sequences where “(a) many of the reads have at least one good, valid alignment, (b) many of the reads are relatively high-quality, and (c) the number of alignments reported per read is small (close to 1)” (“Bowtie”). As a result the Maq algorithm scales better for the longer sequences, but is still 14.9x better than Maq for 76 bp (Langmead et al. 2009). The result will find the location of the alignment.

There is an updated version of Bowtie released in 2012, known as Bowtie 2 (Langmead and Salzberg 2012). Bowtie does not find gapped alignments, however in the recent incarnation of Bowtie, Bowtie 2, it allows for gapped alignments of any number or size of gaps. Bowtie 2 is not a replacement for Bowtie and Bowtie is still superior for bp less than 50. For this project, we will use Bowtie to benchmark against.

## Test Cases to Test Algorithm Limits
For this project, we will use both existing genomes and sequences as well as generate synthetic sequences. Artificial sequences and references sequences will be generated with a known level of error. Error includes both insertions and deletions as well as substitutions that account for between 10-75% of the sequence. 

Both Bowtie and our project’s algorithm will be benchmarked against accuracy. Test data will be generated with a range of errors in the sequences. The synthetic data will be generated with the errors 5%, 10%, 20%, 30%, 40%, 50%, 60%, 70%, and 80%. Bowtie is most effective for sequences less than 75 bp (avg. 50 bp).  We will generate 50 sequences with bp with multiples of 10: 10, 20, …100 bp. From the 50 variation of a sequence of the same length with the same level of error, the average accuracy will be recorded and charted.

This will test the limits of the both algorithms. There will also be sequences of all the lengths with 0% (should allows find) and 100% (should never find) error for edge case testing. 0% and 100% will also be useful for measure a baseline for both algorithms when it can find the sequence the fastest and the slowest. Accuracy is how often the algorithm is able to find the sequence.

As a result, there will be 11 types of errors, with 10 different lengths with 50 variations. Benchmarking will produce a total of 5500 (11*10*50) of sequences that each algorithm will be benchmarked against with the average accuracy stored and charted along with running time. These sequences will be tested against a static reference genome that is 1000 bp long.

Additional artificial sequences will be generated to further test the limits of both algorithms. Additional sequences will includes duplicates of the aligned value to see how both algorithms behave with more than one option, which is a feature that this project’s algorithm can be optimized for.

Sequences being aligned and referenced will be constant between both algorithms to ensure consistent benchmarking between them. In order to test our application as well as Bowtie we will randomly generate DNA sequence data. We will first test our application against the smaller sequences for accuracy and precision. Then we will assess efficiency using the larger sequences from actual DNA sequences retrieved and generated from NCBI.

## References
BEAR. "Sej917/BEAR." GitHub. Sej917, 10 Feb. 2017. Web. 07 Apr. 2017. <https://github.com/sej917/BEAR>.

"Benchmarker for Evaluating the Effectiveness of RNA-Seq Software (BEERS)." Benchmarker for Evaluating the Effectiveness of RNA-Seq Software (BEERS). BEER, n.d. Web. 07 Apr. 2017. <http://www.cbil.upenn.edu/BEERS/>.

"Bowtie." Bowtie: Manual. Sourceforge, n.d. Web. 04 Apr. 2017.

EAGLE. "Sequencing/EAGLE." GitHub. Sequencing, 23 Dec. 2016. Web. 07 Apr. 2017. <https://github.com/sequencing/EAGLE>.

Escalona, Merly, Sara Rocha, and David Posada. "A Comparison of Tools for the Simulation of Genomic Next-generation Sequencing Data." Nature Reviews. Genetics. U.S. National Library of Medicine, 17 Aug. 2016. Web. 05 Apr. 2017.

Langmead, Ben, Cole Trapnell, Mihai Pop, and Steven L. Salzberg. "Ultrafast and Memory-efficient Alignment of Short DNA Sequences to the Human Genome." Genome Biology 10.3 (2009): Web

Langmead, Ben, and Steven L. Salzberg. "Fast Gapped-read Alignment with Bowtie 2." Nature Methods 9.4 (2012): 357-59. Web.

Paolo Ferragina, and Giovanni Manzini. "Opportunistic data structures with applications." Proceedings of the 41st Annual IEEE Symposium on Foundation of Computer Science IEEE (2000): 390-398. Web.
Renaud, G., M. Kircher, U. Stenzel, and J. Kelso. "FreeIbis: An Efficient Basecaller with Calibrated Quality Scores for Illumina Sequencers." Bioinformatics 29.9 (2013): 1208-209. Web.
