---
layout: page
nav:
  - name: Description
    link: "#description"
  - name: Usage
    link: "usage.html"
  - name: Download
    link: "https://github.com/AlgoLab/PIntron/tarball/master"
  - name: Installation
    link: "install.html"
  - name: License
    link: "#license"
  - name: Authors
    link: "#authors"
  - name: Known Issues
    link: "https://github.com/AlgoLab/PIntron/issues"
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

Please cite the following paper when using **PIntron**:
Yuri Pirola, Raffaella Rizzi, Ernesto Picardi, Graziano
Pesole, Gianluca Della Vedova and Paola Bonizzoni.
*PIntron: A fast method for gene structure prediction via maximal
pairings of a pattern and a text,*
BMC Bioinformatics 2012, 13(Suppl 5):S2 [doi:10.1186/1471-2105-13-S5-S2](http://dx.doi.org/10.1186/1471-2105-13-S5-S2)

The program is developed and heavily tested on Linux, but it should
run also on OS X.

## Download

The latest source code version can be downloaded from the [GitHub page](https://github.com/AlgoLab/PIntron) as a [zip](https://github.com/AlgoLab/PIntron/zipball/master) or as [tar](https://github.com/AlgoLab/PIntron/tarball/master) archive.
It is also possible to clone the source repository using the following command:

    git clone https://github.com/AlgoLab/PIntron.git

PIntron is only distributed as source code.

We release only stable versions, therefore you are encouraged to
upgrade always to the latest version.

Anyway, older versions are also [available](https://github.com/AlgoLab/PIntron/tags).

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
