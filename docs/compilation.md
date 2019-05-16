# Introduction [![Build Status](https://travis-ci.org/choishingwan/PRSice.svg?branch=master)](https://travis-ci.org/choishingwan/PRSice)
Here is the guideline for anyone who might want to compile PRSice from source. 

# Prerequisites
For the C++ executable
1. GCC version 4.8.1 or higher
2. CMake version 3.1 or higher (Optional)
3. Git (Optional)
!!! note
    Only the C++ executable need to be built

# Using CMake
With CMake, you can simply do the following:
```
git clone https://github.com/choishingwan/PRSice.git
cd PRSice
mkdir build
cd build
cmake ../
make
```
Then the PRSice executable will be located within PRSice/bin

If you don't have git installed, you can still do
```
curl https://codeload.github.com/choishingwan/PRSice/tar.gz/2.1.3.beta > PRSice.tar.gz
tar -xvf PRSice.tar.gz
cd PRSice-2.1.3.beta
mkdir build
cd build
cmake ../
make
```

!!! Note
    The above procedure was not tested on Windows

# Without CMake
Without CMake, you can simply do the following
```
git clone https://github.com/choishingwan/PRSice.git
cd PRSice
g++ -std=c++11 -O3 -DNDEBUG -march=native -isystem lib -I inc src/*.cpp -lpthread -lz -o PRSice
```
Then PRSice will be located in the current directory

Alternatively, if you don't have git installed, you can still do
```
curl https://codeload.github.com/choishingwan/PRSice/tar.gz/2.2.0 > PRSice.tar.gz
tar -xvf PRSice.tar.gz
cd PRSice-2.2.0.beta
g++ -std=c++11 -O3 -DNDEBUG -march=native -isystem lib -I inc src/*.cpp -lpthread -lz -o PRSice
```

# Intel MKL
If you know how to setup the Intel $\circledR$ MKL library, you can compile PRSice with it to speed up the processing speed. You can use [this](https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/) to help you with the linking.