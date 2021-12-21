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
+ Python 3<
+ Conda 4.7<

## Install - (Mac/Linux/Windows)
### Downloading SeqWho and setting up conda environment
This download example will place SeqWho in your home (~) directory if you are using a linux system. 
Change the ~ to whichever directory you desire if this is not the behavior you want.

```bash
git clone https://github.com/DaehwanKimLab/seqwho ~/seqwho

cd seqwho

conda env create -f base_conda.yml
```

### Adding SeqWho to PATH
Add the above directory (SeqWho) to your PATH environment variable
(e.g. ~/.bashrc) to make the binaries built above and other python scripts
available everywhere:

```bash
$ export PATH=~/seqwho:$PATH
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
> Comma-separated list of files, individual file, or pattern matching syntax to classify. The files can be gzipped or uncompressed.
> For paired-end reads, please treat them as single-end reads.
> Example: `-f file1.fq.gz,file2.fq.gz` or `-f ~/files/*.fq`

### Optional Arguments

* **\-o / \--out** | *Default* : SeqWho_call 
> Directory name for output calls of SeqWho
> Example: `-o SeqWho_call`

