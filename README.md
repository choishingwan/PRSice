# PRSice
PRSice (pronounced 'precise') is a software package for calculating, applying, evaluating and plotting the results of polygenic risk scores (PRS). 
PRSice can run at high-resolution to provide the best-fit PRS as well as provide results calculated at broad P-value thresholds, illustrating results corresponding to either, can thin SNPs according to linkage disequilibrium and P-value ("clumping"), and can be applied across multiple traits in a single run.

Based on a permutation study we estimate a significance threshold of P = 0.001 for high-resolution PRS analyses - the work on this is included in our Bioinformatics paper on PRSice.

PRSice is a software package written in R and C++. PRSice runs as a command-line program with a variety of user-options and is freely available for download below, compatible for Unix/Linux/Mac OS.

For more details on the authors, see: [Jack's homepage](https://kclpure.kcl.ac.uk/portal/en/persons/jack-euesden(972d61b2-89c6-4777-8969-7d88b0c0ece5).html), [Cathryn's homepage](http://www.kcl.ac.uk/lsm/research/divisions/gmm/archive/clusters/bse/lewis/clewis.aspx), [Paul's homepage](http://www.pauloreilly.info/).

## Prerequisite
GCC version 4.8.1 or higher (for c++11)
[Eigen C++](http://eigen.tuxfamily.org/index.php?title=Main_Page) (included)
[Boost library](http://www.boost.org/) (included)

##Installation
You can directly download the binary files for Linux here.
If you want to install PRSice, all you have to do is (The binary file will located in PRSice)
```bash
git clone https://github.com/choishingwan/PRSice.git
cd PRSice
make
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

## MAC users
If you are trying to compile PRSice in MAC without using cmake, then you will need to use
```bash
git clone https://github.com/choishingwan/PRSice.git
cd PRSice
g++ --std=c++11 -I inc/ -isystem lib/ -L /usr/lib/x86_64-redhat-linux5E/lib64 -DNDEBUG -O2 -pthread -Wno-c++11-narrowing src/*.cpp src/*.c -o PRSice
```

### Rosalind users
You can compile a static version using the following command
```bash
git clone https://github.com/choishingwan/PRSice.git
cd PRSice
make
```

### Citation 
If you PRSice in any published work, please cite both the software (as an electronic resource/URL):

Package: PRSice [version]

Authors: Jack Euesden, Cathryn Lewis & Paul O'Reilly

URL: https://github.com/choishingwan/PRSice

and the manuscript(s):

Jack Euesden  Cathryn M. Lewis  Paul F. O’Reilly (2015) PRSice: Polygenic Risk Score software . Bioinformatics 31 (9): 1466-1468

As PRSice is heavily rely on PLINK, please also cite:

Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience, 4.


### Note to Self
PLINK PRS range is inclusive. e.g. 0 - 0.5 includes also SNPs with p-value of 0 and 0.5
