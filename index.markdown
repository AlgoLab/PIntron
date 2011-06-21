---
layout: default
title: PIntron
---

A novel pipeline for gene-structure prediction based on spliced alignment of transcript sequences (ESTs and mRNAs) against a genomic sequence.

## Description

PIntron is composed by four steps: Firstly, alternative alignments of expressed sequences to a reference genomic sequence are implicitly computed and represented in a graph (called embedding graph) by a novel fast spliced alignment procedure. Secondly, biologically meaningful alignment are extracted. Then, a consensus gene structure induced by the previously computed alignments is determined based on a parsimony principle. Finally, the resulting introns are reconciliated and classified according to general biological criteria.

A detailed description of the pipeline is available [here](http://arxiv.org/abs/1005.1514).

Please refer to the file [dist-docs/README.md](https://github.com/AlgoLab/PIntron/blob/master/dist-docs/README.md) for the detailed usage notes.


## Download

The latest stable version of PIntron is **v1.2.0** (*June 21, 2011*).

PIntron is distributed as source code and as pre-built binary packages for some common platforms.

The latest source code version can be downloaded from the [GitHub page](https://github.com/AlgoLab/PIntron) as a [zip](https://github.com/AlgoLab/PIntron/zipball/master) or as [tar](https://github.com/AlgoLab/PIntron/tarball/master) archive.
It is also possible to clone the source repository using the following command:

    git clone https://github.com/AlgoLab/PIntron.git


The binary packages are available from the [AlgoLab webpage](http://www.algolab.eu/PIntron).


## Install

Please refer to the file [dist-docs/INSTALL.md](https://github.com/AlgoLab/PIntron/blob/master/dist-docs/INSTALL.md).


## License

PIntron is distributed under the terms of the GNU Affero General Public License (AGPL), either version 3 of the License, or (at your option) any later version.

PIntron is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.

A copy of the GNU Affero General Public License is distributed along with PIntron and it is also available at the [GNU website](http://www.gnu.org/licenses/).


## Authors

The joint main contributors of this project are Yuri Pirola (<yuri.pirola@gmail.com>) and Raffaella Rizzi (<rizzi@disco.unimib.it>).

The work has also greatly benefited from the contribution and the supervision of Paola Bonizzoni (<bonizzoni@disco.unimib.it>) and Gianluca Della Vedova (<gianluca.dellavedova@unimib.it>).

For any information, please contact one of the authors.


