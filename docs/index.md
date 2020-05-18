<h1>PRSice-2: Polygenic Risk Score software</h1>

PRSice (pronounced 'precise') is a Polygenic Risk Score software for calculating, applying, evaluating and plotting the results of polygenic risk scores (PRS) analyses. Some of the features include:


PRSice (pronounced 'precise') is a Polygenic Risk Score software for calculating, applying, evaluating and plotting the results of polygenic risk scores (PRS) analyses. Some of the features include:

1. High-resolution scoring (PRS calculated across a large number of P-value thresholds)
2. Identify Most predictive PRS
3. Empirical P-values output (not subject to over-fitting)
4. Genotyped (PLINK binary) and imputed (Oxford bgen v1.2) data input
5. Biobank-scale genotyped data can be analysed within hours
6. Incorporation of covariates
7. Application across multiple target traits simultaneously
8. Results plotted in several formats (bar plots, high-res plots, quantile plots)
9. PRSet: function for calculating PRS across user-defined pathways / gene sets

# Executable downloads [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3703335.svg)](https://doi.org/10.5281/zenodo.3703335)[![Coverage Status](https://coveralls.io/repos/github/choishingwan/PRSice/badge.svg?branch=master)](https://coveralls.io/github/choishingwan/PRSice?branch=master)
| Operating System | Link |
| -----------------|:----:|
| Linux 64-bit | [v2.3.0](https://github.com/choishingwan/PRSice/releases/download/2.3.0/PRSice_linux.zip) |
| OS X 64-bit | [v2.3.0](https://github.com/choishingwan/PRSice/releases/download/2.3.0/PRSice_mac.zip) |
| Windows 32-bit | Not available |
| Windows 64-bit | Not available |

!!! Note "Latest Update"
 
    - We now support multi-threaded clumping (separated by chromosome)
    - Genotypes will be stored to memory during clumping (increase memory usage, significantly speed up clumping)
    - Will only generate one .prsice file for all phenotypes
        - .prsice file now has additional column call "Pheno"
    - Introduced `--chr-id` which generate rs id based on user provided formula (see [detail](command_detail.md) for more info)
    - Format of `--base-maf` and `--base-info` are now changed to `<name>:<value>` from `<name>,<value>`
    - Fix a bug related to ambiguous allele dosage flipping when `--keep-ambig` is used
    - Better mismatch handling. For example, if your base file only provide the effective allele A without the non-effective allele information, PRSice will now do dosage flipping if your target file has G/C as effective allele and A /T as an non-effective allele (whereas previous this SNP will be considered as a mismatch)
    - Fix bug in 2.2.13 where PRSice won't output the error message during command parsing stage
    - If user provided the `--stat` information, PRSice will now error out instead of trying to look for BETA or OR in the file. 
    - PRSice should now better recognize if phenotype file contains a header
    - various small bug fix

    update log for previous release can be found [here](update_log.md)



!!! Caution

    We have now fixed window problem. But was unable to access the computer that is used for compilation due to COVID. Will try to compile it when we regain access.


!!! Caution 

    PRSet are currently under open beta - results output are reliable but please report any specific problems to our google group (see Support below)3


# R Packages Requirements

To plot graphs, PRSice requires [R](https://www.r-project.org/) (**version 3.2.3+**) installed.

[Additional steps](extra_steps.md) might be required for Mac and Windows users.

!!! Note "Installing required R packages" 

    PRSice can automatically download all required packages, even without administrative right.
    You can specify the install directory using `--dir`. For example

    ```
        Rscript PRSice.R --dir .
    ```

    will install all required packages under the local directory.

# Quick Start
For Quick start use, please refer to [Quick Start](quick_start.md)

!!! tip "List user options"

    You can also type

    ```
        ./PRSice
    ```

    to view all available parameters unrelated to plotting, or

    ```
        Rscript PRSice.R -h
    ```

    to view all available parameters, including those used for plotting

# Output of Results
You can see the expected output of PRSice [here](step_by_step.md#output-of-results)

# Detailed Guide
You can find a more detailed document explaining the input and output of PRSice in [this page](step_by_step.md)

# Full command line options
You can find all command line options of PRSice under the section *Details of PRSice/PRSet*

## Citation
If you use PRSice, then please cite:

!!! important "Citation"

	Choi SW, and O’Reilly PF. "PRSice-2: Polygenic Risk Score Software for Biobank-Scale Data." GigaScience 8, no. 7 (July 1, 2019). [https://doi.org/10.1093/gigascience/giz082](https://doi.org/10.1093/gigascience/giz082).


## Support
This wiki should contain all the basic instruction for the use of PRSice.
Shall you have any problems, please feel free to start an issue [here](https://github.com/choishingwan/PRSice/issues) or visit our [google group](https://groups.google.com/forum/#!forum/prsice).
You can help us to speed up the debug process by including the log file generated by PRSice.

In addition, you can use the search bar in this webpage to search for specific functions. 

## Authors
For more details on the authors, see:

- [Dr Shing Wan Choi](https:choishingwan.github.io)
- [Dr Paul O'Reilly](http://www.pauloreilly.info/)

PRSice-2 and all new functionalities are coded by:

- [Dr Shing Wan Choi](https://choishingwan.github.io)


## Acknowledgement
PRSice is a software package written in C++ (main) and R (plotting).
The code relies partially on those written in PLINK by [Christopher Chang](https://www.cog-genomics.org/software).
Management of BGEN file is based on BGEN lib written by [Gavin Band](https://bitbucket.org/gavinband/bgen).
We also utilize the [Eigen C++](https://eigen.tuxfamily.org) library, the [gzstream](http://www.cs.unc.edu/Research/compgeom/gzstream/) library. 

