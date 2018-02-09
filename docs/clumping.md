# Introduction
Due to Linkage Disequilibrium (LD) between SNPs,
signals in high-LD regions will be over-represented,
leading to bias in the PRS.
To maximize the signals for PRS analysis, an association
information driven pruning algorithm, i.e. Clumping can be
performed.
That is, when a SNP pair is in LD, only the SNP with the
lowest p-value will be retained, this prevents the more associated
SNP being removed from subsequent analysis, thus preserved the signal.

In PRSice, the PLINK clumping algorithm is implemented where
the r^2^ values computed are based on maximum likelihood haplotype
frequency estimates.
In addition, the LD is calculated in all founder samples (ignoring case control
status).

!!! tip

    If you would like to perform the LD estimates on only the controls, you
    can use `--ld` and `--keep`/`--remove`.

!!! important

    When no reference panel is provided, the target samples will be
    used for LD calculation. However, none of the `--ld-*` parameters
    will take effect. The `--ld-*` parameters will only have an effect
    on the `--ld` file. You will need to resupply the target file to the
    `--ld` parameter for the `--ld-*` parameters to take effect

# Commands
- `--clump-kb`
    The distance for clumping in kb. Default: 250

- `--clump-r2`

    The r^2^ threshold for clumping. Default: 0.1

- `--clump-p`

    The p-value threshold use for clumping. Default: 1.

- `--ld` | `-L`

    LD reference file. Use for estimation of LD during clumping.
    If not provided, will use the post-filtered target genotype
    for LD calculation. Support multiple chromosome input.
    Please see [`--target`](target_file.md#target-file-related-parameters) for more information.
    When you target sample is small (e.g. < 500) and
    external panel of the same population is available (e.g. 1000 genome),
    you might want to use the external reference panel
    to improve the LD estimation for clumping.

- `--ld-geno`

    Filter SNPs based on genotype missingness. Must be a value
    between *0.0* and *1.0*.

- `--ld-info`

    Filter SNPs based on info score. Only used for imputed LD reference.
    The INFO score is calculated as the MaCH imputation r-squared value.
    The pseudo code is represented as follow:

```
    m=Mean of expected genotype
    v=variance of expected genotype
    p=m/2
    p_a = 2p(1-p)
    INFO = v/p_a
```

- `--ld-hard-thres`

    Hard threshold for dosage data. Any call less than
    this will be treated as missing.Default: 0.9

- `--ld-keep`

    File containing the sample(s) to be extracted from
    the LD reference file. First column should be FID and
    the second column should be IID. If `--ignore-fid` is
    set, first column should be IID.
    Mutually exclusive from `--ld-remove`. No effect if `--ld` was not provided

- `--ld-maf`

    Filter SNPs based on minor allele frequency (MAF)

    !!! Note

        When perform MAF filtering on dosage data, we will
        calculate the MAF using the hard-coded genotype

- `--ld-remove`

    File containing the sample(s) to be removed from
    the LD reference file. First column should be FID and
    the second column should be IID. If `--ignore-fid` is
    set, first column should be IID.
    Mutually exclusive from `--ld-keep`

- `--ld-type`

    File type of the LD file. Support bed (binary plink)
    and bgen format. Default: bed\n"

- `--no-clump`

    When set, PRSice will not perform clumping. This is useful when you have
    a pre-clumped list of SNPs and want PRSice to only perform the regression
    analysis.

- `--proxy`

    Proxy threshold for index SNP to be considered
    as part of the region represented by the clumped
    SNP(s). e.g. `--proxy 0.8` means the index SNP will
    represent region of any clumped SNP(s) that has
    r^2^=0.8 even if the index SNP does not physically
    locate within the region
