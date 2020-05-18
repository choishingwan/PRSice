From now on, I will try to archive our update log here. 


# 2020-05-18 (v2.3.0)
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