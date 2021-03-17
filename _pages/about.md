---
layout: page
title: About
permalink: /about/
order: 2
share: false
---

With the vast improvements in sequencing technologies and increased number of protocols, sequencing is finding more applications to answer complex biological problems. Thus, the amount of publicly available sequencing data has tremendously increased in repositories such as SRA, EGA, and ENCODE. With any large online database, there is a critical need to accurately document study metadata, such as the source protocol and organism. In some cases, this metadata may not be systematically verified by the hosting sites and may result a propogation of error through future studies. Furthermore, it is often imperitive to assess FASTQ(A) file and read quality prior to running any time intensive pipeline. SeqWho is a program designed to heuristically assess the quality of sequencing files and reliably classify the organism and protocol type. This is done in an alignment-free algorithm that leverages a Random Forest classifier to learn from native biases in k-mer frequencies and repeat sequence identities between different sequencing technologies and species. We have shown that our method can accurately and rapidly distinguish between human and mouse, nine different sequencing technologies, and both together, 98.32%, 97.86%, and 96.38% of the time in high confidence calls respectively. SeqWho is a powerful method for reliably checking the identity of the sequencing files used in any pipeline and illustrates the program’s ability to leverage k-mer biases.
