## Download

PIntron is only distributed as source code.

The latest source code version can be downloaded from the [GitHub page](https://github.com/AlgoLab/PIntron) as a [zip](https://github.com/AlgoLab/PIntron/zipball/master) or as [tar](https://github.com/AlgoLab/PIntron/tarball/master) archive.
It is also possible to clone the source repository using the following command:

    git clone https://github.com/AlgoLab/PIntron.git

We release only stable versions, therefore you are encouraged to
upgrade always to the latest version.

Anyway, older versions are also [available](https://github.com/AlgoLab/PIntron/tags).

## Requirements

`PIntron` requires the following free software:

- Python v3.0 or newer. Python v3.1 or newer is
  recommended
- Perl (tested on Perl v5.10.1 but it should work also with older
  versions)
- The JSON Perl module (available from CPAN). On most platforms, this
  module can be installed by using the command `cpan JSON` with
  administrator/superuser privileges (e.g., by using the command `sudo
  cpan JSON` on MacOS X)
- A recent version of the standard GNU toolchain (gcc, make).

## Compilation and Installation

`PIntron` is only distributed as source code and must be manually built.

The build process is driven by the GNU make utility and can be performed
by the following invocation:

    make dist

The command will produce a compressed archive `dist/pintron-*.tar.gz`
that can be used for the installation or the execution.

The binary package is composed by the directory `bin`, containing all the
executables needed to run PIntron, and the directory `doc`, containing
the documentation and a simple complete example.
For installing PIntron, you should copy the executables of the `bin/`
directory to a directory of the PATH or to a custom directory.
In the second case, you should specify the custom directory during the
PIntron invocation using the `--bin-dir` program option.



