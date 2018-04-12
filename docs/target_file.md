# Introduction

Currently we support two different formats

### PLINK Binary
A target dataset in PLINK binary format must consist of three files: **.bed**, **.bim**, and a **.fam** file - where bed contains the compressed genotype data, bim contains the SNP information and fam contains the family information. Currently we only support a SNP major PLINK input (default output of the latest PLINK program).

The **.bed** and **.bim** file must have the same prefix.

If you want to use a **.fam** file with a different prefix from the **.bed** and **.bim** file, you can do
```
--target <bim bed prefix>,<fam file>
```

!!! warning

    You MUST ensure the fam file contains the correct number of samples

Missing phenotype data can be coded as NA, or -9 for binary traits and NA for quantitative traits.
!!! Note

    -9 will NOT be considered as missing for quantitative traits

If your binary file is separated into individual chromosomes, then you can place an # in place of the chromosome number.
PRSice will automatically substitute # with 1-22

i.e. If your files are chr1.<bed|bim|fam>,chr2.<bed|bim|fam>,...,chr22.<bed|bim|fam>, just use
```
--target chr#
```

!!! Note

    If an external fam file is provided, the substitution will not be performed on it. This is because all chromosome should have the same fam file. fam file from any of the chromosome should be ok

!!! todo

    We will incorporate a new function in next version where we allow different number of autosomes to be specified.
    This should allow users to use PRSice on other organisms.
    However, due to complication in handling non-diploid chromosomes, we will only
    support diploid organisms (and exclude all sex chromosomes), i.e. Mouse

### BGEN
We also support BGEN v1.2. To specify a BGEN file, you simply add `--type bgen` or `--ld-type bgen` to your PRSice command

As BGEN does not store the phenotype information and sometime not even the sample ID, you **must** provide
a phenotype file (`--pheno-file`) for PRSice to run. Or you can provide a sample file using

```
--target <bgen prefix>,<sample file>
```

!!! Note

    We will still require this sample file even if you are doing `--no-regress` as we do want the sample ID.  we might loosen this requirement later on.

With BGEN, a number of other PRSice options become effective:
- `--hard`: Normally, with BGEN format, we will calculate the PRS using the dosage information.
But you can ask for hard-thresholding by using the `--hard` option. Then we will code the SNP as the genotype with
probability larger than threshold specified by `--hard-thres`. If no such genotype is presented, the SNP will be
coded as missing
- `--hard-thres`: The genotype probability threshold. SNPs with no genotype having a probability larger than this
threshold will be treated as missing


!!! Note
    Similar to PLINK files, you can use the # substitution for BGEN files that are separated into different chromosomes

## Phenotype files
PRSice also support an external phenotype file as in input using the `--pheno-file`
This file should contain the FID and IID (or just IID if `--ignore-fid` is set)
and the target phenotype.

You can specify the target phenotype name using `--pheno-col`.
Multiple phenotype can be provided using a comma separated list. e.g.

```
  --pheno-col A,B,C
```

will cause PRSice to perform the analysis on phenotype A, B and C.

!!! Note
    When multiple phenotype is provided, PRSice will generate a set of
    results for each phenotype with the prefix of `<Prefix>.<Pheno>` where
    `<Prefix>` is the parameter provided to the `--out` option

# Target File Related Commands

- `--binary-target`

    Indicate whether the target phenotype is binary or not.
    Either T or F should be provided where T represent a binary phenotype.
    For multiple phenotypes, the input should be separated by comma without space.

    Default: **T** if `--beta` is set and **F** otherwise

- `--geno`

    Filter SNPs based on gentype missingness. Must be a value
    between *0.0* and *1.0*.

- `--info`
    Filter SNPs based on info score. Only used for imputed target.
    The INFO score is calculated as the MaCH imputation r-squared value.
    The pseudo code is represented as follow:

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

        When perform MAF filtering on dosage data, we will
        calculate the MAF using the hard-coded genotype

- `--nonfounders`

    By default, PRSice will exclude all non-founders from the analysis.
    When this flag is set, non-founders will be included in the
    regression model but will still be excluded from LD estimation.

- `--pheno-col`

    Headers of phenotypes to be included from the phenotype file.
    When multiple phenotypes are provided, the phenotype name will be
    used as part of the file output prefix

- `--pheno-file` | `-f`

    Phenotype file containing the phenotype(s).
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
    You can also specify a separate fam/sample file by
    `--target <prefix>,<fam/sample file>`


- `--type`

    File type of the target file. Support bed (binary plink) and bgen format. Default: bed

# Dosage Related Commands
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
