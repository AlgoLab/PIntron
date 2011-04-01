#  PIntron

A novel pipeline for gene-structure prediction based on
spliced alignment of transcript sequences (ESTs and mRNAs) against a genomic sequence.

## Description

This software implements a novel algorithm to efficiently predict the
(basic) gene structure by aligning expressed sequence tags (ESTs) and mRNAs against a genomic sequence.

## Installing PIntron

See `INSTALL.md` file.

## Usage

Running the *pintron* program with suitable options achieves the desired
effect. Running pintron with option *-h* gives a detailed list of all available
options.
The two main options *-g* and *-s* allow to specify a genomic sequence file and a file
of transcripts (ESTs and mRNAs) respectively.

### Input file formats

The genomic sequence file is a FASTA file whose header should adhere to the
following structure

      >chrZZ:start:end:strand

where:

- `ZZ` specifies the chromosome
- `start` and `end` specify the starting and ending positions of the genomic sequence
on the reference chromosome.
- `strand` specifies the strand of the genomic sequence and it must be +1 (5'3' direction)
or -1 (3'5' direction).

The provided nucleotide sequence must be on the direction specified by *strand* in the header
FASTA.
The genomic sequence file may have a FASTA header whose format is
different from the above standard *but*, in this case, rare prediction
errors might arise.

The file of transcripts is a multiFASTA file containing all the transcript sequences (ESTs and
mRNAs). The header of each sequence must contain the following field:

    /gb=XXXXXX

and it should contain the following field:

    /clone_end=YY

where `XXXXXX` is the GenBank accession number (or in any case a unique identifier within the input
  transcript file) and `YY` is the strand which the transcript has been read from. `YY` can
  be either `3'` (for direction 5'3', or strand +1)
  or `5'` (for direction 3'5', or strand -1).
  By default, when the `/clone_end` field is not specified, the transcript is
  considered on strand +1. Furthermore, a transcript
  sequence, that is also a RefSeq sequence, is forced on strand +1.


### List of command-line parameters

PIntron can be run with the following options:

* `-g`	`--genomic=FILE`		FILE containing the genomic sequence (*default=*genomic.txt)
* `-s`	`--est=FILE`            FILE containing the transcript sequences (*default=*ests.txt)
* `-o`	`--output=FILE`	output JSON FILE with ... (*default=*pintron-full-output.json)
* `-t`  `--gtf=FILE` output GTF FILE with the isoforms that have a CDS annotation (*default=*pintron-cds-annotated-isoforms.gtf)
* `-l`	`--logfile=FILE`	output FILE containing a log (*default=*pintron-log.txt)
* `-z`	`--compress`	compress output (*default=*no)
* `-c`	`--do-not-continue`	do not resume a previously interrupted computation (*default=*no)
* `-b`	`--bin-dir=DIRECTORY`	DIRECTORY containing the programs (*default=*system PATH)
* `-n`	`--organism=NAME`	NAME of the organism originating the transcript sequences (*default=*unknown)
* `-e`	`--gene=NAME`	NAME of the gene originating the transcript sequences (*default=*unknown)
* `-k`	`--keep-intermediate-files`	keep all intermediate or temporary files  (*default=*no)
* `--extended-gtf=FILE`	output GTF FILE with all the predicted isoforms (*default=*pintron-all-isoforms.gtf)
* `--set-max-factorization-time=INT` set a time limit (in mins) for the spliced alignment step (*default=*60)
* `--set-max-factorization-memory=INT` set a limit (in MiB) for the memory used by the spliced alignment step (*default=*3000 MiB, approx. 3GB)
* `--set-max-exon-agreement-time=INT` Set a time limit (in mins) for the exon agreement step (*default=*15)
* `--set-max-intron-agreement-time=INT` Set a time limit (in mins) for the intron agreement step (*default=*30)

A typical invocation of PIntron is:

    pintron --genomic=genomic.txt --est=ests.txt -z --output=pintron-full-output.json --logfile=pintron-log.txt

or, equivalently:

    pintron -g genomic.txt -s ests.txt -z -o pintron-full-output.json -l pintron-log.txt



#### Advanced customization for experienced users

The user can choose the parameters to be used in the spliced alignment step, by providing an
optional config file (that must be named *config.ini*) in the current working directory.
The option *-k* allows keeping the dump file *config-dump.ini* of all the parameters actually used
in the spliced alignment step. Both files are composed of rows in the following format:

        <parameter>="<value>"

where:

- `<parameter>` is the name of a parameter
- `<value>` is the provided value.

An example of row is:

        min-factor-length="15"

Some significant parameters that the user can provide are:

* `min-factor-length` The length (nt) of the minimum common transcript-genome factor (*default=*15)
* `min-intron-length` The minimum length (nt) allowed for an intron (*default=*60)
* `max-intron-length`  The maximum length (nt) allowed for an intron (*default=*0, this parameter is ignored)
* `max-prefix-discarded` The maximum length (nt) of a transcript prefix that can be discarded. (*default=*50)
* `max-suffix-discarded` The maximum length (nt) of a transcript suffix that can be discarded. (*default=*50)
* `retain-externals` If *false*, discard the first (and the last, if a polyA chain has not been found) factors of
a spliced alignment. (*default=*true)

### Output

The output file:

* specified with the `--output` option is a JSON file, containing informations
  on the gene structure, the isoforms and the exons. Since JSON is mainly a
  key, value hash, we have chosen keys that are as self-explaining as
  possible. Each exon has keys  "sequence", "relative start", "relative end",
  "chromosome start", "chromosome end" (absolute coordinates), "3utr length",
  "5utr length". Each isoform has keys "canonical start codon?",
  "canonical end codon?", "NMD flag", "polyA?", "annotated CDS?", "start",
  "end", "number exons", "Type", "length", "coding length", "reference?", "from RefSeq?",
  "exons" (a list of exons). Each gene has keys  "length_genomic_sequence",
  "version",  "number_isoforms", "isoforms" (an hash of isoforms with a
  progressive ID has the only key).
* specified with the `--gtf` option is a standard GTF file, containing the isoforms that have a CDS annotation.
* specified with the `--extended-gtf` option is a GTF file, containing all the predicted isoforms.


### A complete example

The *doc/example* directory contains the input and the output of a complete execution of
the PIntron pipeline over the human ENCODE SLC5A1 gene (chromosome 22, plus strand).

The input files are:

* `genomic.txt` file of the ENCODE region ENm004 (nucleotide sequence provided on plus strand)
* `ests.txt` file of the transcript sequences (ESTs/mRNAs) related to the considered gene

The command (from `dist/pintron-latest` if PIntron has been built from
source, or from the directory on which the binary package has been
extracted) for obtaining the output files contained in `doc/example` is:

        ./bin/pintron \
                --bin-dir=bin \
                --genomic=doc/example/genomic.txt \
                --est=doc/example/ests.txt \
                --organism=human \
                --gene=SLC5A1 \
                --output=doc/example/pintron-full-output.json \
                --gtf=doc/example/pintron-cds-annotated-isoforms.gtf \
                --extended-gtf=doc/example/pintron-all-isoforms.gtf \
                --logfile=doc/example/pintron-log.txt
