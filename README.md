## DeepMSA: multiple sequence alignment generation for protein structure prediction ##

#### Installing the package ####

This package is developed for 64bit Linux only. If you use 64bit Linux, you
should be able to use this package without compiling anything.

In theory, the package should also work for Mac OS, provided that you compile
the respective components appropriately. In particular, this package contains
source code and/or binary executables from the following packages.

[1] hhsuite 2.0.16
    Source code and binaries included. Check
    [hhsuite-readme.txt](hhsuite-readme.txt) for compilation instruction.

[2] legacy BLAST 2.2.26
    Binaries included. Check 
    [ncbi-blast-legacy](https://github.com/kad-ecoli/ncbi-blast-legacy)
    for source code and 
    ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/
    for binary executable.

[3] kClust 1.0
    Binaries included. Check [kClust](https://github.com/soedinglab/kClust)
    for source code.

[4] Clustal Omega 1.2.4
    Binary included. Check [CLustal Omega](http://www.clustal.org/omega)
    for source code.

[5] HMMER 3.1b2
    Modified binaries included. Check
    [qhmmer](https://github.com/kad-ecoli/qhmmer) to compile our modified
    source code for the qhmmsearch, qhmmbuild, and qjackhmmer binaries.

[6] MSAParser, our own program to calculate Nf. Check 
    https://github.com/kad-ecoli/MSAParser/ for source code.

#### Running the package ####

The main program is [scripts/build_MSA.py](scripts/build_MSA.py). You can
run the script without any command line argument to get the help message.
<b>The program needs python2.7</b>; python3 will result in error.

#### License ####

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

If you find this program useful, please cite:

C Zhang, W Zheng, SM Mortuza, L Yang, Y Zhang (2019)
DeepMSA: constructing deep multiple sequence alignment to improve contact 
prediction and fold-recognition for distant-homology proteins.
Bioinformatics. DOI:
[10.1093/bioinformatics/btz863](https://doi.org/10.1093/bioinformatics/btz863)
