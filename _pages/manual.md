---
layout: page
title: Manual
permalink: /manual/
order: 3
share: false
---

{: .no_toc}

- TOC
{:toc}

# Introduction

SeqWho has two major functions: 1) model training 2) file classification. In both cases SeqWho function by k-merizing the reads in FASTQ(A) files and counting the 1-5mers to build a frequency array. SeqWho can then utilize either a full HISAT2 repeat index during training or the internal repeat index in the SeqWho.ix file during file classification. SeqWho uses these core functions to train a Random Forest classifier as the core mode of classification. Theoretically, there is no limit to the number of classifications possible in SeqWho. It depends mostly on getting enough data to properly train the models.

---

# SeqWho Set-up

## Requirements
SeqWho is written in Python 3 and will not function in Python 2. While it is certainly possible to manually set-up Python 3 to use SeqWho, we recommend using a conda environment built from the environment.yml included with SeqWho for optimal performance. This is the method we will describe here. 

### Basic Recommended Requirements
+ Python >3
+ Conda >4.7

## Install - (Mac/Linux/Windows)
SeqWho provides two installation method.
 - PIP
 - Downloading SeqWho from GitHub

### PIP
We recommend installing SeqWho in a separate conda environment.

```bash
conda create -n seqwho_v1 python=3.7 pip
conda activate seqwho_v1
pip install seqwho
```

### Downloading SeqWho and setting up conda environment
This download example will place SeqWho in your home (~) directory if you are using a linux system. 
Change the ~ to whichever directory you desire if this is not the behavior you want.

```bash
git clone https://github.com/DaehwanKimLab/seqwho ~/seqwho

cd seqwho

conda env create -f base_conda.yml
```

#### Adding SeqWho to PATH
Add the above directory (SeqWho) to your PATH environment variable
(e.g. ~/.bashrc) to make the binaries built above and other python scripts
available everywhere:

```bash
$ export PATH=~/seqwho/seqwho_lib:$PATH
```

---

# Running SeqWho 

We have separated this section into model building and file classification

Please activate the conda environment before running SeqWho

```bash
$ conda activate seqwho_v1
```

## Model Building
Usage:
```bash
$ seqwho_buildindex.py -r [REPEAT FILES] -l [LABEL FILE] [OPTIONS]
```

*Note: this building step requires all training files and supporting files be in the current working directory where seqwho_buildindex is being run*

### Required Arguments

* **\-r / \--repeats** | *Default* : *None* 
> A comma-separated list of repeat indicies to use. Typically from HISAT2 
> Example: `-r mouse.rep,human.rep`

* **\-l / \--labels** | *Default* : *none* 
> CSV file with labels and file names for training. Format: file_name,species,sequence_type  
> Example: `-l labels.csv`

### Optional Arguments

* **\-m / \--mask** | *Default* : (empty) 
> A comma-separated list of any file types to omit in the training
> Example: `-m whole_genome_sequencing,rnaseq`

* **\-k / \--ksize** | *Default* : 5
> Set the max size of the k-mers to count
> Example: `-k 5` represents 1-5 mers

* **\-j / \--repeat-kmer** | *Default* : 31  
> Size of repeat k-mers to use. Default is 31-mer  
> Example: `-j 31`

* **\-o / \--out** | *Default* : SeqWho.ix  
> Name of output index file
> Example: `-o SeqWho`

* **\--rebuild** | *Default* : `FALSE` 
> Overwrite any existing index if one is detected with the same name as determined in -o

* **\-v / \--verbose** | *Default* : `FALSE`
> Show output messages

### Generate Repeat Files for seqwho_buildindex.py
To generate the repeat files for seqwho_buildindex.py, you need to install [HISAT2] and passing the genome file to `hisat2-repeat`.  
Example:
```bash
$ hisat2-repeat mouse_genome.fa mouse
```
The file `mouse.rep.100-65535.seed` is the repeat file for seqwho_buildindex.py `-r` option.

[HISAT2]:https://daehwankimlab.github.io/hisat2/main/
## File Calling
Usage:
```bash
$ seqwho.py -x [SEQWHO INDEX] -f [FILE(S) ... ] [OPTIONS]
```

### Required Arguments

* **\-x / \--index** | *Default* : *None* 
> Path to SeqWho index file and name of file
> Example: `-x SeqWho.ix`

* **\-f / \--files** | *Default* : *none*
> Space-separated list of files, individual file, or pattern matching syntax to classify. The files can be gzipped or uncompressed.
> For paired-end reads, please treat them as single-end reads.
> Example: `-f file1.fq.gz file2.fq.gz` or `-f ~/files/*.fq`

### Optional Arguments

* **\-o / \--out** | *Default* : SeqWho_call 
> Directory name for output calls of SeqWho
> Example: `-o SeqWho_call`

# Results and output files
SeqWho generates three kinds of files:

1) PNG containing four plots, one with single nucleotide quality core counts, the second with averag quality score by position, the third with cound of reads with a givent length, and the final a heat map of nucleotide frequency by position in the read.

2) A Json file with all of the data neede to generate the plots in 1, and the following data:

 - Estimated Read Number: Number of read estimated to be in the file
 - File format: File format (fasta or fastq)
 - Biased 5' end 6-mers": Any 6-mers at the 5' end that are more frequent than expected by chance
 - Mean Quality: Average quality of bases
 - Mean Read Len: Average length of reads
 - Perc PolyA: Percent reads that are poly-A
 - Perc GC Cont: Percent GC content
 - Perc Reads w N: Percent of reads with an N
 - Reads Passed: Number of reads passing into model
 - Reads Omitted: Number of reads filtered because of low quality or are poly-A
 - Perc Passed: Percent of reads that passes to model testing
 - Call: Call of the file


3) TSV with all of the information in 2 above

A PNG is generated for each file tested and dropped into a separate results folder  while only one JSON and TSV file are generated with all of the results imbedded.