---
layout: page
toc: false
---

## Description

**PIntron** is a novel pipeline for gene-structure prediction based on spliced alignment of transcript sequences (ESTs and mRNAs) against a genomic sequence.

**PIntron** is composed by four steps:
Firstly, alternative alignments of expressed sequences to a reference
genomic sequence are implicitly computed and represented in a graph
(called *embedding graph*) by a novel fast spliced alignment procedure.
Secondly, biologically meaningful alignment are extracted.
Then, a consensus gene structure induced by the previously computed
alignments is determined based on a parsimony principle.
Finally, the resulting introns are reconciliated and classified
according to general biological criteria.

The program is developed and heavily tested on Linux, but it should
run also on OS X.

## Features

## Screenshots
**PIntron** is a command-line application that can be used on local
workstations or on remote servers.
Anyway, the computed results can be used with third-party applications
to visualize and further process the results.

Here, for example, is the set of full-length isoforms (upper part)
computed by **PIntron** for gene
[TP53](http://www.ncbi.nlm.nih.gov/gene/7157) and visualized by the
[UCSC Genome Browser](http://genome.ucsc.edu/).

[![img](images/pintron-v1.2.28-TP53.png)](images/pintron-v1.2.28-TP53.png)

