From now on, I will try to archive our update log here. 

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
