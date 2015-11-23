# Maligner

Maligner is a tool for aligning molecular or insilico restriction maps to a reference map. Maligner comes with several different modes of alignment:

 - `maligner_dp`: Uses dynamic program and allows global-local alignments of a query against a reference. Allows for unmatched sites in both the query and reference.

 - `malign_ix` : Uses a more restrictive but faster mode of index based alignment.

 - `malign_vd` : Allows for partial prefix or suffix alignments of a query against a reference, which can be used to find split alignments.
 

## Installation

Maligner is built using [cmake](https://cmake.org/download/), which can be installed using your package manager. Building Maligner requires a C++ compiler C++11 support.

To build with cmake,  

```bash
git clone https://github.com/LeeMendelowitz/maligner.git maligner
cd maligner # Make sure you are in the repo directory
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
make install
```

This will install compiled binaries `maligner_dp`, `maligner_ix`, and `maligner_vd` and additional python utility scripts into the directory `build/bin`.

The `malignpy` python package is installed to 'build/lib'. Many of the Maligner utility scripts for working with maps files and alignment files depend on `malignpy`. In order to use these scripts, you must symlink malignpy into your working directory or modify your `PYTHONPATH` environment variable:

```bash
export PYTHONPATH=/path/to/build/lib:$PYTHONPATH
```

### Dependencies

Building Maligner requires a C++ compiler with C++11 support. The build has been
tested with g++ 4.8.3 on Red Hat Linux and Apple LLVM 6.0 on Mac OS X.

Maligner also comes with several python scripts for working
with and converting alignment files. These scripts require the following python libraries: numpy, scipy, pandas, and BioPython.

You can install these dependencies using [pip](https://pip.pypa.io/en/stable/):

```bash
pip install -r requirements.txt
```

## Getting Started

 - See the [ecoli example script](https://github.com/LeeMendelowitz/maligner/blob/master/examples/ecoli_example.sh) in the repository.

 - See the [wiki](https://github.com/LeeMendelowitz/maligner/wiki)

 - Ask for command line help: `maligner_dp --help`

