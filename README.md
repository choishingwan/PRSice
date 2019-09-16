# PRSice [![Build Status](https://travis-ci.org/choishingwan/PRSice.svg?branch=master)](https://travis-ci.org/choishingwan/PRSice)
PRSice (pronounced 'precise') is a software package for calculating, applying, evaluating and plotting the results of polygenic risk scores (PRS). 
PRSice can run at high-resolution to provide the best-fit PRS as well as provide results calculated at broad P-value thresholds, illustrating results corresponding to either, can thin SNPs according to linkage disequilibrium and P-value ("clumping"), and can be applied across multiple traits in a single run.

Based on a permutation study we estimate a significance threshold of P = 0.001 for high-resolution PRS analyses - the work on this is included in our Bioinformatics paper on PRSice.

PRSice is a software package written in R and C++. PRSice runs as a command-line program with a variety of user-options and is freely available for download below, compatible for Unix/Linux/Mac OS

## NOTE
Please refer to our [website](https://www.prsice.info/) for more update instructions

## Prerequisite
GCC version 4.8.1 or higher (for c++11)
R version 3.2.3 or higher (for plotting)

## Installation
You can directly download the binary files [here](https://github.com/choishingwan/PRSice/releases/).
If you want to install PRSice, all you have to do is (The binary file will located in PRSice)
```bash
git clone https://github.com/choishingwan/PRSice.git
cd PRSice
g++ --std=c++11 -I inc/ -isystem lib/ -DNDEBUG -O3 -march=native src/*.cpp -lz -lpthread -o PRSice
```
Or if you have CMake version 3.1 or higher, you can do (The binary file will located in PRSice/bin)
```bash
git clone https://github.com/choishingwan/PRSice.git
cd PRSice
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
make
```

### Citation 
If you PRSice in any published work, please cite the following manuscript:

Choi SW, and O’Reilly PF. "PRSice-2: Polygenic Risk Score Software for Biobank-Scale Data." GigaScience 8, no. 7 (July 1, 2019). https://doi.org/10.1093/gigascience/giz082.

### Note to Self
PLINK PRS range is inclusive. e.g. 0 - 0.5 includes also SNPs with p-value of 0 and 0.5
