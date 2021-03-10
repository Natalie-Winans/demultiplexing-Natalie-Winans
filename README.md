# Demultiplexing

This project was completed as part of the Bioinformatics and Genomics Master's Program at the University of Oregon.

The objective was to design an algorithm to quality filter and demultiplex short-read sequencing data and assess levels of index hopping. Demultiplexing is required when sequencing libraries are pooled -- the resulting reads must be sorted by their barcodes or indexes. Index hopping can interfere with accurate downstream analyses by causing incorrect or ambiguous read classification.

The program `Winans-Demultiplexer.py` can be called from the command line and requires the following arguments:

```
-r1, --read1        FASTQ file containing read 1 ("forward reads")
-r2, --read2        FASTQ file containing read 2 ("index 1")
-r3, --read3        FASTQ file containing read 3 ("index 2")
-r4, --read4        FASTQ file containing read 4 ("reverse reads")
-i, --index_file    Tab-separated file containing index data
```

The index file must be in the following format:
```
sample	group	treatment	index	index sequence
1	2A	control	B1	GTAGCGTA
2	2B	control	A5	CGATCGAT
3	2B	control	C1	GATCAAGG
4	2C	mbnl	B9	AACAGCGA
6	2D	mbnl	C9	TAGCCATG
7	2E	fox	C3	CGGTAATC
8	2F	fox	B3	CTCTGGAT
10	2G	both	C4	TACCGGAT
11	2H	both	A11	CTAGCTCA
14	3B	control	C7	CACTTCAC
15	3C	mbnl	B2	GCTACTCT
16	3D	mbnl	A1	ACGATCAG
17	3E	fox	B7	TATGGCAC
19	3F	fox	A3	TGTTCCGT
21	3G	both	B4	GTCCTAAG
22	3H	both	A12	TCGACAAG
23	4A	control	C10	TCTTCGAC
24	4A	control	A2	ATCATGCG
27	4C	mbnl	C2	ATCGTGGT
28	4D	mbnl	A10	TCGAGAGT
29	4E	fox	B8	TCGGATTC
31	4F	fox	A7	GATCTTGC
32	4G	both	B10	AGAGTCCA
34	4H	both	A8	AGGATAGC
```

The program output consists of six FASTQ files and a summary text file.

The output FASTQ files include:
1. Forward reads with matched indexes 
2. Reverse reads with matched indexes
3. Forward reads with hopped indexes
4. Reverse reads with hopped indexes
5. Forward reads with unknown or low quality indexes
6. Reverse reads with unknown or low quality indexes

Example summary file:
```
--------------------------
| Demultiplexing Summary |
--------------------------

Number of matched reads: 304980270
Percent matched reads: 83.96%
Number of hopped reads: 517612
Percent hopped reads: 0.14%
Number of unknown reads: 57748853
Percent unknown reads: 15.9%
Total reads: 363246735


Matched counts per sample:
--------------------------
SAMPLE	INDEX SEQ	COUNT	PERCENT OF MATCHED READS
01	GTAGCGTA	7450201	  2.44%
02	CGATCGAT	5225776	  1.71%
03	GATCAAGG	6085915	  2.0%
04	AACAGCGA	8178191	  2.68%
06	TAGCCATG	9852258	  3.23%
07	CGGTAATC	4498136	  1.47%
08	CTCTGGAT	32163349    10.55%
10	TACCGGAT	69307073	22.73%
11	CTAGCTCA	16162895	5.3%
14	CACTTCAC	3833640	  1.26%
15	GCTACTCT	6610857	  2.17%
16	ACGATCAG	7441721	  2.44%
17	TATGGCAC	10195805	3.34%
19	TGTTCCGT	14786868	4.85%
21	GTCCTAAG	8164223	  2.68%
22	TCGACAAG	3548541	  1.16%
23	TCTTCGAC	39149148	12.84%
24	ATCATGCG	9264615	  3.04%
27	ATCGTGGT	6357656	  2.08%
28	TCGAGAGT	10658212	3.49%
29	TCGGATTC	4163314	  1.37%
31	GATCTTGC	3425453	  1.12%
32	AGAGTCCA	10378366	3.4%
34	AGGATAGC	8078057	  2.65%


Counts per hopped index pair:
-----------------------------
INDEX PAIR	COUNT	PERCENT OF HOPPED READS
('GTAGCGTA', 'CGATCGAT')	144	0.03%
('GTAGCGTA', 'GATCAAGG')	156	0.03%
('GTAGCGTA', 'AACAGCGA')	222	0.04%
('GTAGCGTA', 'TAGCCATG')	197	0.04%
('GTAGCGTA', 'CGGTAATC')	104	0.02%

.                           .   .
.                           .   .
.                           .   .
```


