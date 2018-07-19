# Introduction
Here, we detail some of the decisions we made during the implementatino of PRSice

# Support of BGEN v1.3
We have purposefully disabled support to BGEN v1.3 to avoid the inclusion of the zstd library. 
This is because
- UKBB is v1.2
- We are not familiar with the licensing of zstd library (developed by facebook)

# Removal of PCA calculation
The main goal of PRSice 2 is to support the polygenic score analysis on large scale data. 
With such data, the calculation of PCA on the fly will be time consuming and will require specific algorithms
such as those implemented in [flashPCA](https://github.com/gabraham/flashpca).
In order to support the in-place PCA calculation, not only will we have to implement the flashPCA algorithm, we 
will also need to implement prunning, which is required prior to PCA calculation. 

Due to the lack of man power, we therefore decided that we will not implement the PCA calculation. 
Another reasoning is that we believe users should first examine the PCA results before directly applying them
to the PRS analysis. 

