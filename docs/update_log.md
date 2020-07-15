From now on, I will try to archive our update log here. 
# 2020-07-15 ([v2.3.2](https://github.com/choishingwan/PRSice/tree/e4b146e7d118277660fdfc3f5813eaebe61433ce))
- Fix off by one error in PRSet best score output
- Fix problem for bgen file when sample selection is performed on bgen files containing sample information


# 2020-05-30 ([v2.3.1.e](https://github.com/choishingwan/PRSice/tree/dc04caa38c05e9f15484edaeb1cfce341da7cf1d))
- Fix bug where SNPs without missingness will be wrongly considered as having 100% missingness
- Fix error log where PRSice should now correct stat if a parameter is missing the required arguments 

# 2020-05-29 ([v2.3.1.d](https://github.com/choishingwan/PRSice/tree/1201ef4e3811dbe099fc9d49b7a463f48dc6025c))
- Fix segmentation fault when `--ld` is used

# 2020-05-28 ([v2.3.1.c](https://github.com/choishingwan/PRSice/tree/56b84ea6051cba23ee91bbbad4ebde272582bbd6))
- Fix problem with missing covariate
- Fix Rscript such that it properly read in phenotype file when `--pheno-co
l` is specified

# 2020-05-26 ([v2.3.1.b](https://github.com/choishingwan/PRSice/tree/9e756f9c8fe9f9ed24c4b5e6c770f64f3112eeb1))
- Fix best score output when `--ignore-fid` is used
- Also fix Rscript covariate and phenotype file read when handling IDs star
t with 00 and when `--ignore-fid` is used

# 2020-05-26 ([v2.3.1.a](https://github.com/choishingwan/PRSice/tree/86b002170316c63d1b7255c5f1d5f136242802c0))
- Fix bar plot with covariate. Was plotting the full R2 instead of the PRS.R2

# 2020-05-23 ([v2.3.1](https://github.com/choishingwan/PRSice/tree/91f4265ad5c30643c0676c6bb37a404fff021bc3))
- Update Rscript such that it match features in executable (thus avoid problem in plotting)
- Fix a bug where PRSice will crash when there are missing covariates


# 2020-05-21 ([v2.3.0.e](https://github.com/choishingwan/PRSice/tree/2b057f0eafa28762ec0c1245bc2f20aacadda05b))
- Fix Rscript bar plot problem

# 2020-05-21 ([v2.3.0.d](https://github.com/choishingwan/PRSice/tree/8784ab58b5171c5e4bbc5341de5baa68f5f5238f))
- Fix problem introduced by previous fix.
- Was hoping 2.3.0's unit test will help reducing the amount of bugs. Sorry for the troubles.

# 2020-05-20 ([v2.3.0.c](https://github.com/choishingwan/PRSice/tree/3fca49456ea5f0d84e01c06d0c491fbb5917181a))
- Fix all score output format
- Fix problem with `--no-regress`. Might still have problem with `--no-regress --score con-std`
    

# 2020-05-19 ([v2.3.0.b](https://github.com/choishingwan/PRSice/tree/a999a862b83599497bcea3fa16cde340dca52e11))
- Fix error where sample selection will distort phenotype loading, loading the wrong phenotype to wrong sample. As
this is a major bug, we deleted the previous 2 releases. Sorry for the troubles.

# 2020-05-19 ([v2.3.0.a](https://github.com/choishingwan/PRSice/tree/87c8571f8b27d39cfe6d8ec3b00e059d0ecf0376))
- Fix output error where we always say 0 valid phenotype were included for continuous trait
- Fix problem with permutation where PRSice will crash if input are rank deficient
- Fix problem when provide a binary phenotype file with a fam file containing -9 as phenotype, PRSice will wrongly state that there are no phenotype presented
- Fix problem in Rscript where if sample ID is numeric and starts with 0, the best file will not merge with the phenotype file, causing 0 valid PRS to be observed

# 2020-05-18 ([v2.3.0](https://github.com/choishingwan/PRSice/tree/2.3.0))
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
