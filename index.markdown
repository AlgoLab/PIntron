---
layout: page
toc: false
---

## Description

**PIntron** is a novel pipeline for gene-structure prediction based on spliced alignment of transcript sequences (ESTs and mRNAs) against a genomic sequence.

The program is developed and heavily tested on Linux, but it should
run also on OS X.

## Features

Starting from a genomic sequence and the ESTs sequence, encoded as
MultiFASTA file, we produce:

*  the exons and introns
*  the full-length isoforms


## Screenshots
**PIntron** is a command-line application that can be used on local
workstations or on remote servers.
Anyway, the computed results can be used with third-party applications
to visualize and further process the results.

Here, for example, is the set of full-length isoforms (upper part)
computed by **PIntron** for gene
[TP53](http://www.ncbi.nlm.nih.gov/gene/7157) and visualized by the
[UCSC Genome Browser](http://genome.ucsc.edu/).

[![img](images/pintron-v1.2.28-TP53-hs437460.png)](images/pintron-v1.2.28-TP53-hs437460.png)

## Main ideas

**PIntron** is composed of steps:
1. Alternative alignments of expressed sequences to a reference
genomic sequence are computed and represented in a graph
(called *embedding graph*), using an ad-hoc fast spliced alignment
procedure.
2. Biologically meaningful alignment are filtered.
3. A draft consensus gene structure is constructed
4. Introns are reconciliated and classified 
according to general biological criteria.
