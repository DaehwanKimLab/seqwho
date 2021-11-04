# SeqWho - An accruate and rapid FASTQ(A) file origin classifier

This is the official SeqWho Repository.

SeqWho is a reliable and extremely rapid program designed to determine a FASTQ(A) sequencing file identity, both source protocol and species of origin. This is accomplished using an alignment-free algorithm that leverages a Random Forest classifier that learns from biases in k-mer frequencies and repeat sequence identity. SeqWho is capable of achieving greater than 96% accuracy in its ability to classify files.

You can find the Documentation for SeqWho at:
https://daehwankimlab.github.io/seqwho/

## First time setup
SeqWho is written in Python 3 and we recommend using a conda environment built from the environment.yml included with SeqWho for optimal performance.

Please read https://daehwankimlab.github.io/seqwho/manual/ for more details.

### Download pre-trained SeqWho index

| Species | Libraries | Index |
| --- | --- | --- |
| Human, Mouse | Amplicon, ChIP-Seq, WGS, WES, miRNA-Seq, RNA-Seq, Bisulfite-Seq, DNase-Seq, ATAC-Seq | [Index]( https://cloud.biohpc.swmed.edu/index.php/s/sP48taKmymSkJBM/download)  <sub>[md5sum](https://cloud.biohpc.swmed.edu/index.php/s/9bk57S65LycK5ts/download)</sub> | 
| Human, Mouse, Rattus norvegicus | ChIP-Seq, WGS, RNA-Seq | [Index](https://cloud.biohpc.swmed.edu/index.php/s/B6dLibRmjqkWbNY/download) <sub>[md5sum](https://cloud.biohpc.swmed.edu/index.php/s/GcCyaFqeaQSCsTd/download)</sub> |
| Human, Mouse, Arabidopsis thaliana, Caenorhabditis elegans, Drosophila melanogaster, Saccharomyces cerevisiae | ChIP-Seq, WGS, RNA-Seq | [Index](https://cloud.biohpc.swmed.edu/index.php/s/ys6Qa87cY2HyJEJ/download) <sub>[md5sum](https://cloud.biohpc.swmed.edu/index.php/s/7AZdAHEc6iRBYSP/download)</sub> |


## Current release
v1.0.3 - Added option to select number of reads drawn from files during model building

v1.0.2 - Removed extra commas in some fields to facilitate CSV conversion

v1.0.1 - Addition of test files and scripts

v1.0.0 - Initial public release

