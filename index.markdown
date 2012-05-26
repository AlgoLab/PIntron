---
layout: default
title: PIntron
---

A novel pipeline for gene-structure prediction based on spliced alignment of transcript sequences (ESTs and mRNAs) against a genomic sequence.

## Description

**PIntron** is composed by four steps:
Firstly, alternative alignments of expressed sequences to a reference
genomic sequence are implicitly computed and represented in a graph
(called *embedding graph*) by a novel fast spliced alignment procedure.
Secondly, biologically meaningful alignment are extracted.
Then, a consensus gene structure induced by the previously computed
alignments is determined based on a parsimony principle.
Finally, the resulting introns are reconciliated and classified
according to general biological criteria.

A detailed description of the pipeline is available
[here](http://arxiv.org/abs/1005.1514).

Please cite the following paper when using **PIntron**:   
Yuri Pirola, Raffaella Rizzi, Ernesto Picardi, Graziano
Pesole, Gianluca Della Vedova and Paola Bonizzoni.
*PIntron: A fast method for gene structure prediction via maximal
pairings of a pattern and a text,* 
BMC Bioinformatics 2012, 13(Suppl 5):S2 [doi:10.1186/1471-2105-13-S5-S2](http://dx.doi.org/10.1186/1471-2105-13-S5-S2)

Please refer to this [page](usage.html) or to the
[`README`](https://github.com/AlgoLab/PIntron/blob/master/dist-docs/README.md)
file in the source distribution for the detailed usage notes.


## Download

The latest stable version of PIntron is **v1.2.17** (*Nov. 15, 2011*).

PIntron is only distributed as source code.

The latest source code version can be downloaded from the [GitHub page](https://github.com/AlgoLab/PIntron) as a [zip](https://github.com/AlgoLab/PIntron/zipball/master) or as [tar](https://github.com/AlgoLab/PIntron/tarball/master) archive.
It is also possible to clone the source repository using the following command:

    git clone https://github.com/AlgoLab/PIntron.git

Older versions are also [available](https://github.com/AlgoLab/PIntron/tags).


## Install

Please refer to this [page](install.html) or to the
[`INSTALL`](https://github.com/AlgoLab/PIntron/blob/master/dist-docs/INSTALL.md)
file in the source distribution for the detailed installation instructions.


## License

PIntron is distributed under the terms of the GNU Affero General Public
License (AGPL), either version 3 of the License, or (at your option) any
later version.

PIntron is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public
License for more details.

A copy of the GNU Affero General Public License is distributed along
with PIntron and it is also available at the
[GNU website](http://www.gnu.org/licenses/).


## Authors

The joint main contributors of this project are Yuri Pirola
(<yuri.pirola@gmail.com>) and Raffaella Rizzi (<rizzi@disco.unimib.it>).

## Known Issues

We use the issue tracker of
[GitHub](https://github.com/AlgoLab/PIntron/issues) to report all
software bugs as well as the progress on their solutions.

Please report [there](https://github.com/AlgoLab/PIntron/issues) any problem you find with PIntron.
