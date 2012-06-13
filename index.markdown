---
layout: page
toc: false
---

## Description

**PIntron** is a novel pipeline for gene-structure prediction based on
spliced alignment of transcript sequences (ESTs and mRNAs) against a
genomic sequence.

The program is developed and heavily tested on Linux, but it should
run also on OS X.

## Features

Starting from a genomic sequence and the transcript (EST and/or mRNA)
sequences, encoded as (Multi)FASTA files, we predict:

*  the exons and introns of the gene(s)
*  the full-length isoforms

The output can be easily inspected by the users and can be also easily
parsed by simple scripts written in any programming languages.
Moreover, the predictions are also reported in the standard GTF format
and, thus, can be easily used and visualized by other commonly used
software systems (such as the [UCSC Genome Browser][]).


## Screenshots
**PIntron** is a command-line application that can be used on local
workstations or on remote servers.
Anyway, the computed results can be used with third-party applications
to visualize and further process the results.

Here, for example, is the set of full-length isoforms (upper part)
computed by **PIntron** for gene
[TP53](http://www.ncbi.nlm.nih.gov/gene/7157) and visualized by the
[UCSC Genome Browser][].

[![img](images/pintron-v1.2.28-TP53-hs437460.png)](images/pintron-v1.2.28-TP53-hs437460.png)

## Main ideas

**PIntron** is a pipeline composed of five steps:

1. Alternative alignments of expressed sequences to a reference genomic
   sequence are computed and represented in a graph (called *embedding
   graph*), using a new ad-hoc fast spliced alignment procedure.
2. Biologically meaningful alignments are filtered.
3. A draft consensus gene structure is constructed.
4. Introns are reconciliated and classified according to general
   biological criteria.
5. A minimal set of full-length isoforms that explains the reconstructed
   introns and the spliced alignments is finally computed using the
   method described [here](http://dx.doi.org/10.1089%2fcmb.2008.0028).



[UCSC Genome Browser]: http://genome.ucsc.edu/
