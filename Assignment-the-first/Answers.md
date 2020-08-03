# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read 1 |
| 1294_S1_L008_R2_001.fastq.gz | index 1 |
| 1294_S1_L008_R3_001.fastq.gz | index 2 |
| 1294_S1_L008_R4_001.fastq.gz | read 2 |

2. Per-base NT distribution
    1. For Read 1 ("forward"):\
    ![R1 Mean Scores](https://github.com/2020-bgmp/demultiplexing-Natalie-Winans/blob/master/Assignment-the-first/1294_S1_L008_R1_001.fastq_mean_scores.png?raw=true)
    2. For Read 2 (Index 1):\
    ![R2 Mean Scores](https://github.com/2020-bgmp/demultiplexing-Natalie-Winans/blob/master/Assignment-the-first/1294_S1_L008_R2_001.fastq_mean_scores.png?raw=true)
    3. For Read 3 (Index 2):\
    ![R3 Mean Scores](https://github.com/2020-bgmp/demultiplexing-Natalie-Winans/blob/master/Assignment-the-first/1294_S1_L008_R3_001.fastq_mean_scores.png?raw=true)
    4. For Read 4 ("reverse")"\
    ![R4 Mean Scores](https://github.com/2020-bgmp/demultiplexing-Natalie-Winans/blob/master/Assignment-the-first/1294_S1_L008_R4_001.fastq_mean_scores.png?raw=true)
    
## Part 2
1. Define the problem\
Need to take reads from 4 input FASTQ files (read 1, read 2, index 1, and index2), append index pairs to the headers of "forward" and "reverse" reads, and sort reads into three categories based on indices: matched, hopped, and unknown. 
2. Describe output
* Within the matched category, there will be 48 files: for 24 indices, one file each for the "forward" and "reverse" reads.
* There will be 2 output files each for the hopped and unkown categories.
* This gives a total of 52 output files with sorted reads.
* Output should also include stats: number of matched indices (per index pair and total), number of hopped indices, and number of unknown indices. This can be provided in a summary output file.
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [4 expected output FASTQ files](../TEST-output_FASTQ).
4. [Pseudocode](https://github.com/2020-bgmp/demultiplexing-Natalie-Winans/blob/master/demux_pseudocode.txt)
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
