#!/usr/bin/env python

import argparse
import gzip
import itertools as it

def get_args():
        parser = argparse.ArgumentParser(description="A program to find the mean quality score by base position from a fastq file")
        parser.add_argument("-r1", "--read1", help='file containing read 1 ("forward")')
        parser.add_argument("-r2", "--read2", help="file containing read 2 (index 1)")
        parser.add_argument("-r3", "--read3", help="file containing read 3 (index 2)")
        parser.add_argument("-r4", "--read4", help='file containing read 4 ("reverse")')
        parser.add_argument("-i", "--index_file", help="name of file containing index data")
        return parser.parse_args() 
args = get_args()
R1 = args.read1
R2 = args.read2
R3 = args.read3
R4 = args.read4 
index_file = args.index_file

# R1 = '/Users/nataliewinans/bioinformatics/Bi622/Demux/demultiplexing-Natalie-Winans/TEST-input_FASTQ/test_R1_in.fq.gz'
# R2 = '/Users/nataliewinans/bioinformatics/Bi622/Demux/demultiplexing-Natalie-Winans/TEST-input_FASTQ/test_R2_in.fq.gz'
# R3 = '/Users/nataliewinans/bioinformatics/Bi622/Demux/demultiplexing-Natalie-Winans/TEST-input_FASTQ/test_R3_in.fq.gz'
# R4 = '/Users/nataliewinans/bioinformatics/Bi622/Demux/demultiplexing-Natalie-Winans/TEST-input_FASTQ/test_R4_in.fq.gz'
# index_file = '/Users/nataliewinans/bioinformatics/Bi622/Demux/demultiplexing-Natalie-Winans/indexes.txt'

def reverse_comp(string):
    """Takes DNA sequence (string) and returns reverse complement"""
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join(comp[i] for i in reversed(string))

def convert_phred(x):
    """Converts ASCII code to decimal quality score"""
    return ord(x) - 33

def mean_score(line):
    """Takes string of ASCII scores and returns numerical mean"""
    sum = 0
    for char in line:
        score = convert_phred(char)
        sum += score
    return(sum/len(line))

def get_indices(index_file):
    """Creates dictionary with keys = index sequences (e.g. 'TCGAGAGT') and values = sample number (e.g. '3')"""
    with open(index_file, 'rt') as fh:
        index_dict = {}
        for i,line in enumerate(fh):
            if i > 0:
                seq = [x.rstrip() for x in line.split('\t')][4] #list comprehension - removes \r and \t
                sample = line.split('\t')[0]
                index_dict[seq] = sample.zfill(2) #adds leading zero to single digits for later sorting
    return(index_dict)

def get_index_seq(sample):
    """Takes the sample number and returns the index sequence"""
    for seq, val in index_dict.items():
        if val == sample:
            return seq

index_dict = get_indices(index_file)
#print(index_dict)

#create dictionary to hold index pair permutations and counts
poss_pairs = list(it.permutations(index_dict.keys(), 2))
pdict = {}
for pair in poss_pairs:
    pdict[pair] = 0

#initialize record counters for each output file
#dictionary to hold matched counts
matched_counts = {}
for sample in index_dict.values():
    matched_counts[sample] = 0

unknown_count = 0
hopped_count = 0
tot_matched_count = 0
total_reads = 0

#initialize lists to hold one record from each input file
R1_rec = []
R2_rec = []
R3_rec = []
R4_rec = []

#open output files
file_dict_R1 = {}
file_dict_R4 = {}
for sample in index_dict.values():
    file_dict_R1[sample] = open('sample_%s_matched_R1.fastq' % sample, 'w')
    file_dict_R4[sample] = open('sample_%s_matched_R4.fastq' % sample, 'w')

hopped_R1 = open('hopped_R1_out.fastq', 'w')
hopped_R4 = open('hopped_R4_out.fastq', 'w')
unknown_R1 = open('unknown_R1_out.fastq', 'w')
unknown_R4 = open('unknown_R4_out.fastq', 'w')

#open and loop over 4 input files simultaneously
with gzip.open(R1,'rt') as R1, gzip.open(R2, 'rt') as R2, gzip.open(R3, 'rt') as R3, gzip.open(R4, 'rt') as R4:
    n=4
    for line_R1, line_R2, line_R3, line_R4 in zip(R1, R2, R3, R4):
        #populate record lists
        R1_rec.append(line_R1.rstrip())
        R2_rec.append(line_R2.rstrip())
        R3_rec.append(line_R3.rstrip())
        R4_rec.append(line_R4.rstrip())
        if len(R1_rec) == n:
            #add indices from R2 and R3 to header lines in corresponding records in R1 and R4
            R1_rec[0] = R1_rec[0] + ' ' + R2_rec[1] + '-' + R3_rec[1]
            R4_rec[0] = R4_rec[0] + ' ' + R2_rec[1] + '-' + R3_rec[1]
            total_reads += 1
            #compare indices and write records to output files
            if 'N' in R2_rec[1] \
                or 'N' in R3_rec[1] \
                or mean_score(R2_rec[3]) < 30 \
                or mean_score(R3_rec[3]) < 30 \
                or R2_rec[1] not in index_dict \
                or reverse_comp(R3_rec[1]) not in index_dict:
                unknown_R1.write('\n'.join(R1_rec) + '\n')
                unknown_R4.write('\n'.join(R4_rec) + '\n')
                unknown_count += 1
            
            elif R2_rec[1] == reverse_comp(R3_rec[1]): #matched records
                sample = index_dict[R2_rec[1]]
                file_dict_R1[sample].write('\n'.join(R1_rec) + '\n')
                file_dict_R4[sample].write('\n'.join(R4_rec) + '\n')
                tot_matched_count += 1
                if sample in matched_counts:
                    matched_counts[sample] += 1

            else:
                hopped_R1.write('\n'.join(R1_rec) + '\n')
                hopped_R4.write('\n'.join(R4_rec) + '\n')
                hopped_count += 1
                #increment count for index pair permutation
                if (R2_rec[1], reverse_comp(R3_rec[1])) in pdict: 
                    pdict[R2_rec[1], reverse_comp(R3_rec[1])] += 1

            #empty record lists
            R1_rec = []
            R2_rec = []
            R3_rec = []
            R4_rec = []

#close all output files
for index in index_dict.values():
    file_dict_R1[index].close()
    file_dict_R4[index].close()
hopped_R1.close()
hopped_R4.close()
unknown_R1.close()
unknown_R4.close()

#write summary stats to output file
with open('demux_summary.txt', 'w') as fh:
    fh.write("--------------------------\n")
    fh.write("| Demultiplexing Summary |\n")
    fh.write("--------------------------\n\n")
    fh.write("Number of matched reads: " + str(tot_matched_count) + '\n')
    fh.write("Percent matched reads: " + str(round(float(tot_matched_count)/float(total_reads)*100,2)) + '%\n')
    fh.write("Number of hopped reads: " + str(hopped_count) + '\n')
    fh.write("Percent hopped reads: " + str(round(float(hopped_count)/float(total_reads)*100,2)) + '%\n')
    fh.write("Number of unknown reads: " + str(unknown_count) + '\n')
    fh.write("Percent unknown reads: " + str(round(float(unknown_count)/float(total_reads)*100,2)) + '%\n')
    fh.write("Total reads: " + str(total_reads) + '\n')
    fh.write("\n\nMatched counts per sample:\n")
    fh.write("--------------------------\n")
    fh.write("SAMPLE\tINDEX SEQ\tCOUNT\tPERCENT OF MATCHED READS\n")
    for sample in sorted(matched_counts):
        fh.write(str(sample) + '\t' + get_index_seq(sample) + '\t' + str(matched_counts[sample]) + '\t' + str(round(matched_counts[sample]/tot_matched_count*100, 2)) + '%' + '\n')
    fh.write("\n\nCounts per hopped index pair:\n")
    fh.write("-----------------------------\n")
    fh.write("INDEX PAIR\tCOUNT\tPERCENT OF HOPPED READS\n")
    for pair, count in pdict.items():
        fh.write(str(pair) + '\t' + str(count) + '\t' + str(round(count/hopped_count*100, 2)) + '%' '\n')






