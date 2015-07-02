# Maligner

Maligner is a tool for aligning molecular or insilico restriction maps to a reference map. Maligner comes with several different modes of alignment:

 - `maligner_dp`: Uses dynamic program and allows global-local alignments of a query against a reference. Allows for unmatched sites in both the query and reference.

 - `malign_ix` : Uses a more restrictive but faster mode of index based alignment.

 - `malign_sv` : Allows for partial prefix or suffix alignments of a query against a reference, which can be used to find split alignments.
 

## Installation

Maligner can be built using cmake. The build requires a C++ compiler C++11 support.
From the repository directory:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
make install
```

## Dependencies

Building Maligner requires a C++ compiler with C++11 support. The build has been
tested with g++ 4.8.3 on Red Hat Linux and Apple LLVM 6.0 on Mac OS X.

Maligner also comes with several Python scripts for working
with and converting alignments files. These scripts require the following python libraries: numpy, pandas, lxml, xml, and BioPython. 






