#  PIntron

A novel pipeline for gene-structure prediction based on
spliced alignment of transcript sequences (ESTs and mRNAs) against a genomic sequence.


## Requirements

`PIntron` requires the following free software:

- Python v3.0 or newer (available from `http://www.python.org/download/`)
- Perl (tested on Perl v5.10.1 but it should work also with older
  versions)
- The JSON Perl module (available from CPAN). On most platforms, this
  module can be installed by using the command `cpan JSON` with
  administrator/superuser privileges (e.g., by using the command `sudo
  cpan JSON` on MacOS X)


## Compilation and Installation

`PIntron` is only distributed as source code and must be manually built.

The build process is driven by the GNU make utility and can be performed
by the following invocation:

    make dist

For the compilation process, you will need the standard build tools such
as the aforementioned GNU make utility and a recent C compiler (tested on
`gcc v4.4`).

The command will produce a compressed archive `dist/pintron-*.tar.gz`
that can be used for the installation or the execution.

The binary package is composed by the directory `bin`, containing all the
executables needed to run PIntron, and the directory `doc`, containing
the documentation and a simple complete example.
For installing PIntron, you should copy the executables of the `bin/`
directory to a directory of the PATH or to a custom directory.
In the second case, you should specify the custom directory during the
PIntron invocation using the `--bin-dir` program option.



