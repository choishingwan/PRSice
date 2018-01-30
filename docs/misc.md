# Introduction

Below are some other parameters available for PRSice

# Command
- `--all-score`
    Output PRS for ALL threshold.
    !!! warning
        This will generate a huge file

- `--exclude`

    File contains SNPs to be excluded from the analysis.
    Mutually exclusive from `--extract`

- `--extract`

    File contains SNPs to be included in the analysis.
    Mutually exclusive from `--exclude`

- `--ignore-fid`

    Ignore FID for all input. When this is set,
    first column of all file will be assume to
    be IID instead of FID

- `--logit-perm`

    When performing permutation, still use logistic
    regression instead of linear regression. This
    will substantially slow down PRSice.

    !!! note

        One problem with using `--logit-perm` is that some of
        the permuted phenotype might be suffer from perfect
        separation. This leads to the GLM logistic model not
        being able to be converge (thus terminating PRSice).

        If you encounter such problem, you might want to exclude
        the `--logit-perm` option. In most case, the p-value of the
        linear model should be similar to the logistic model

- `--keep-ambig`

    Keep ambiguous SNPs. Only use this option
    if you are certain that the base and target
    has the same A1 and A2 alleles

- `--out` | `-o`

    Prefix for all file output.

    !!! Note

        If multiple target phenotypes are included (e.g. using `--pheno-col`),
        the phenotype will be appended to the output prefix

        If multiple gene set are included, the name of the set will be appended
        to the output prefix (after the phenotype (if any))

- `--perm`

    Number of permutation to perform. This will
    generate the empirical p-value. Recommend to
    use value larger than 10,000

    !!! note

        Essentially, when permutation is required, PRSice will perform the following
        operation (with slight difference in implementation for optimum performance)

        1. Perform normal PRSice across all thresholds and obtain p-value of the
        most significant threshold
        2. Repeat PRSice analysis *N* times with permuted phenotype. Count the
        number of time where the p-value of the most significant threshold for
        the permuted

- `--seed` | `-s`

    Seed used for permutation. If not provided,
    system time will be used as seed. This will
    allow the same results to be generated when
    the same seed and input is used

- `--print-snp`

    Print all SNPs used to construct the best PRS

- `--thread` | `-n`

    Number of thread use

    !!! tip

        You can ask for maximum number of thread by using `--thread max`

    !!! note

        Multi-threading will actually be slower if you require more thread than
        what is available. To avoid this, PRSice will limit the maximum number
        of thread used to the number of core available on the system as detected by PRSice.

- `--help` | `-h`

    Display the help messages
