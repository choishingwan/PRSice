# PRSice [![Build Status](https://travis-ci.org/choishingwan/PRSice.svg?branch=master)](https://travis-ci.org/choishingwan/PRSice)
PRSice (pronounced 'precise') is a software package for calculating, applying, evaluating and plotting the results of polygenic risk scores (PRS). 
PRSice can run at high-resolution to provide the best-fit PRS as well as provide results calculated at broad P-value thresholds, illustrating results corresponding to either, can thin SNPs according to linkage disequilibrium and P-value ("clumping"), and can be applied across multiple traits in a single run.

Based on a permutation study we estimate a significance threshold of P = 0.001 for high-resolution PRS analyses - the work on this is included in our Bioinformatics paper on PRSice.

PRSice is a software package written in R and C++. PRSice runs as a command-line program with a variety of user-options and is freely available for download below, compatible for Unix/Linux/Mac OS.

For more details on the authors, see: [Jack's homepage](https://kclpure.kcl.ac.uk/portal/en/persons/jack-euesden(972d61b2-89c6-4777-8969-7d88b0c0ece5).html), [Cathryn's homepage](http://www.kcl.ac.uk/lsm/research/divisions/gmm/archive/clusters/bse/lewis/clewis.aspx), [Paul's homepage](http://www.pauloreilly.info/).

## NOTE
Please refer to our [WIKI](https://github.com/choishingwan/PRSice/wiki) for more update instructions

## Prerequisite
GCC version 4.8.1 or higher (for c++11)
R version 3.2.3 or higher (for plotting)

## Installation
You can directly download the binary files for Linux [here](https://github.com/choishingwan/PRSice/releases/).
If you want to install PRSice, all you have to do is (The binary file will located in PRSice)
```bash
git clone https://github.com/choishingwan/PRSice.git
cd PRSice
git checkout beta_testing
g++ --std=c++11 -I inc/ -isystem lib/ -lz -DNDEBUG -O2 -pthread src/*.cpp -o PRSice
```
Or if you have CMake version 3.1 or higher, you can do (The binary file will located in PRSice/bin)
```bash
git clone https://github.com/choishingwan/PRSice.git
cd PRSice
git checkout beta_testing
mkdir build
cd build
cmake ../
make
```

### Rosalind users
You can compile a static version using the following command
```bash
git clone https://github.com/choishingwan/PRSice.git
cd PRSice
git checkout beta_testing
make
```

### Citation 
If you PRSice in any published work, please cite both the software (as an electronic resource/URL):

Package: PRSice [version]

Authors: Shing Wan Choi, Jack Euesden, Cathryn Lewis & Paul O'Reilly

URL: https://github.com/choishingwan/PRSice

and the manuscript(s):

Jack Euesden  Cathryn M. Lewis  Paul F. O’Reilly (2015) PRSice: Polygenic Risk Score software . Bioinformatics 31 (9): 1466-1468

### Note to Self
PLINK PRS range is inclusive. e.g. 0 - 0.5 includes also SNPs with p-value of 0 and 0.5
