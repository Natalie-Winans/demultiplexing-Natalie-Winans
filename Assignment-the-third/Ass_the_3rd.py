#!/usr/bin/env python

import matplotlib.pyplot as plt
import argparse
#from pathlib import Path
import gzip
import itertools as it

def get_args():
        parser = argparse.ArgumentParser(description="A program to find the mean quality score by base position from a fastq file")
        parser.add_argument("-f", "--filename", help="name of input file")
        parser.add_argument("-i", "--index_file", help="name of file containing indexes")
        parser.add_argument("-l", "--seq_len", help="length of quality score line", type=int)
        return parser.parse_args() 
args = get_args()
#file = args.filename
#index_file = args.index_file

R1 = '/Users/nataliewinans/bioinformatics/Bi622/Demux/demultiplexing-Natalie-Winans/TEST-input_FASTQ/test_R1_in.fq'
R2 = '/Users/nataliewinans/bioinformatics/Bi622/Demux/demultiplexing-Natalie-Winans/TEST-input_FASTQ/test_R2_in.fq'
R3 = '/Users/nataliewinans/bioinformatics/Bi622/Demux/demultiplexing-Natalie-Winans/TEST-input_FASTQ/test_R3_in.fq'
R4 = '/Users/nataliewinans/bioinformatics/Bi622/Demux/demultiplexing-Natalie-Winans/TEST-input_FASTQ/test_R4_in.fq'
index_file = '/Users/nataliewinans/bioinformatics/Bi622/Demux/demultiplexing-Natalie-Winans/indexes.txt'

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
    """Create dictionary with keys = index sequences and values = index codes"""
    with open(index_file, 'rt') as fh:
        index_dict = {}
        for i,line in enumerate(fh):
            if i > 0:
                seq = [x.rstrip() for x in line.split('\t')][4] #list comprehension - removes \r and \t
                index_dict[seq] = line.split('\t')[3]
    return(index_dict)

index_dict = get_indices(index_file)
#print(index_dict)

#create dictionary to hold index pair permutations and counts
poss_pairs = list(it.permutations(index_dict.keys(), 2))
pdict = {}
for pair in poss_pairs:
    pdict[pair] = 0

#open output files
file_dict_R1 = {}
file_dict_R4 = {}
for index in index_dict.values():
    file_dict_R1[index] = open('%s_matching_R1.fastq' % index, 'w')
    file_dict_R4[index] = open('%s_matching_R4.fastq' % index, 'w')

hopped_R1 = open('hopped_R1_out.fastq', 'w')
hopped_R4 = open('hopped_R4_out.fastq', 'w')
unknown_R1 = open('unknown_R1_out.fastq', 'w')
unknown_R4 = open('unknown_R4_out.fastq', 'w')

#initialize lists to hold one record from each input file
R1_rec = []
R2_rec = []
R3_rec = []
R4_rec = []

#initialize record counters for each output file
h_r1 = 0
h_r4 = 0
u_r1 = 0
u_r4 = 0

#open and loop over 4 input files simultaneously
with open(R1,'r') as R1, open(R2, 'r') as R2, open(R3, 'r') as R3, open(R4, 'r') as R4:
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
            #compare indices and write records to output files
            if 'N' in R2_rec[1] \
                or 'N' in R3_rec[1] \
                or mean_score(R2_rec[3]) < 20 \
                or mean_score(R3_rec[3]) < 20 \
                or R2_rec[1] not in index_dict \
                or reverse_comp(R3_rec[1]) not in index_dict:
                unknown_R1.write('\n'.join(R1_rec))
                unknown_R4.write('\n'.join(R4_rec))
                u_r1 += 1
                u_r4 += 1
            elif R2_rec[1] == reverse_comp(R3_rec[1]):
                #write R1_rec to matched_R1
                #write R4_rec to matched_R4
                print(R1_rec)
                print(R4_rec)
            else:
                hopped_R1.write('\n'.join(R1_rec))
                hopped_R4.write('\n'.join(R4_rec))
                h_r1 += 1
                h_r4 += 1
                if (R2_rec[1], reverse_comp(R3_rec[1])) in pdict: #increments count for index pair permutation
                    pdict[R2_rec[1], reverse_comp(R3_rec[1])] += 1
                    

                
         
    

            #empty record lists
            R1_rec = []
            R2_rec = []
            R3_rec = []
            R4_rec = []


print(h_r1)
print(h_r4)
print(u_r1)
print(u_r4)




#close output files
# for index in index_dict.values():
#     file_dict_R1[index] = open('%s_matching_R1.fastq' % index, 'w')
#     file_dict_R4[index] = open('%s_matching_R4.fastq' % index, 'w')
hopped_R1.close()
hopped_R4.close()
unknown_R1.close()
unknown_R4.close()


