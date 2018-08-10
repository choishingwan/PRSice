# Basic Commands
- `--binary-target`

    Indicate whether the target phenotype is binary or not.
    Either T or F should be provided where T represent a binary phenotype.
    For multiple phenotypes, the input should be separated by comma without space.

    Default: **F** if `--beta` is set and **T** otherwise

- `--geno`

    Filter SNPs based on gentype missingness. Must be a value
    between *0.0* and *1.0*.

- `--info`
    Filter SNPs based on info score. Only used for imputed target data.
    The INFO score is calculated as the MaCH imputation r-squared value, 
    represented by the following pseudo code:
```
    m=Mean of expected genotype
    v=variance of expected genotype
    p=m/2
    p_a = 2p(1-p)
    INFO = v/p_a
```

- `--keep`

    File containing the sample(s) to be extracted from
    the target file. First column should be FID and
    the second column should be IID. If `--ignore-fid` is
    set, first column should be IID
    Mutually exclusive from `--remove`

 - `--maf`

    Filter SNPs based on minor allele frequency (MAF)

    !!! Note

        When perform MAF filtering on dosage data, the MAF 
        is calculated using the hard-coded genotype

- `--nonfounders`

    By default, PRSice will exclude all non-founders from the analysis.
    When this flag is set, non-founders will be included in the
    regression model but will still be excluded from LD estimation.

- `--pheno-col`

    Headers of phenotypes to be included from the phenotype file.
    When multiple phenotypes are provided, the phenotype name will be
    used as part of the file output prefix

- `--pheno-file` | `-f`

    Tab or space delimited phenotype file containing the phenotype(s).
    First column must be FID of the samples and
    the second column must be IID of the samples.
    When `--ignore-fid` is set, first column must
    be the IID of the samples. Must contain a
    header if `--pheno-col` is specified


- `--prevalence` | `-k`

    Prevalence of all binary trait.
    If provided, PRSice will adjust the ascertainment bias of the R2.

    !!! Note
        When multiple binary trait is found,
        prevalence information must be provided for
        all of them (Either adjust all binary traits,
        or don't adjust at all).
        For example, if there are 3 traits A, B and C,
         where A and C are binary traits with population
         prevalence of 0.1 and 0.2 respectively. The correct input should be
         `--binary-target T,F,T --prevalence 0.1,0.2`

- `--remove`

    File containing the sample(s) to be removed from
    the target file. First column should be FID and
    the second column should be IID. If `--ignore-fid` is
    set, first column should be IID
    Mutually exclusive from `--keep`

- `--target` | `-t`

    Target genotype file. Currently support
    both BGEN and binary PLINK format. For
    multiple chromosome input, simply substitute
    the chromosome number with #.
    PRSice will automatically replace # with 1-22.
    A separate fam/sample file can be specified by
    `--target <prefix>,<fam/sample file>`


- `--type`

    File type of the target file. Support bed (binary plink) and bgen format. Default: bed

# Dosage Related Commands
- `--allow-inter`

    Allow the generate of intermediate file. This will
    speed up PRSice when using dosage data as clumping
    reference and for hard coding PRS calculation

- `--hard-thres`

    Hard threshold for the dosage data. SNPs with be coded as
    the genotype with probability higher than the threshold
    specified here. If no such genotype is presented, the
    SNP will be treated as missing.

    !!! Note

        If dosage data is used as a LD reference, it will always be
        hard coded to calculate the LD

        Default: 0.9

- `--hard`

    When set, will use hard thresholding instead of dosage for PRS construction.
    Default is to use dosage.
