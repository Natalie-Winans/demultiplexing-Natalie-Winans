#!/usr/bin/env python

import matplotlib.pyplot as plt
import argparse
from pathlib import Path
import gzip


def get_args():
        parser = argparse.ArgumentParser(description="A program to find the mean quality score by base position from a fastq file")
        parser.add_argument("-f", "--filename", help="name of input file")
        parser.add_argument("-l", "--seq_len", help="length of quality score line", type=int)
        return parser.parse_args() 
args = get_args()
file = args.filename
seq_len = args.seq_len

scores = []    
for i in range(seq_len):
    scores.append(0.0)

with gzip.open(file) as fh:
    LN = 0
    for line in fh:
        LN += 1
        line = line.strip('\n')
        if LN % 4 == 0:
            for pos, char in enumerate(line):
                qual_score = ord(char) - 33
                scores[pos] += qual_score #sums quality scores for all reads at each base position
NR = LN/4
for pos, sum in enumerate(scores):
    mean = sum/NR
    scores[pos] = mean
    #list 'scores' now holds mean score at each base position

file_name = Path(file).stem

fig=plt.figure(dpi= 100)
plt.bar(range(seq_len), scores, align='center', color='green')
plt.xlabel("Base Position")
plt.ylabel("Mean Quality Score")
plt.title("Mean Quality Score by Base Pair Index")
plt.savefig(file_name + '_mean_scores.png')



