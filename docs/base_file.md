- `--A1`

    Column header containing the **effective allele**.
    There isn't any standardized label for the effective allele, therefore
    extra care must be taken to ensure the correct label is provided, otherwise, 
    the effect will be flipped. 

- `--A2`

    Column header containing **non-effective allele**.

- `--base` | `-b`

    Base (i.e. GWAS) association file. This is a whitespace delimited file
    containing association results for SNPs on the base phenotype.
    This file can be *gzipped* (must have the .gz suffix).
    For PRSice to run, the base file must contain the effective allele
    (`--A1`), effect size estimates (`--stat`), p-value for association
    (`--pvalue`), and the SNP ID (`--snp`).

- `--beta`

    This flag is used to indicate if the test statistic is in the form
    of BETA. When not set, PRSice assume the statistic is in the form
    of Odd Ratios, perform natural log transformation on the test statistic.

- `--bp`

    Column header containing the coordinate of SNPs. When provided,
    the coordinate of the SNPs will be scrutinized between the base and target file.
    SNPs with mismatched coordinate will be excluded.

- `--chr`

    Column header containing the chromosome information. When provided,
    the chromosome information of the SNPs will be scrutinized between the base
    and target file. SNPs with mismatched chromosome information will be automatically excluded.

- `--index`

    If set, assume the base columns are INDEX instead of the name of the corresponding
    columns. Index should be 0-based (start counting from 0)

- `--info-base`

    Base INFO score filtering. Format should be `<Column name>,<Threshold>`.
    SNPs with info score less than `<Threshold>` will be ignored.
    It is useful to perform INFO score filtering to remove SNPs
    with low imputation confidence score. By default, PRSice will search for
    the **INFO** column in your base file and perform info score filtering with
    threshold of **0.9**. You can disable this behaviour by using `--no-default`

- `--maf-base`
    Base minor allele frequency (MAF) filtering.
    Format should be `<Column name>,<Threshold>`.
    SNPs with MAF less than `<Threshold>` will be ignored.
    Additional column can be provided (e.g. different filtering threshold
    for case and control), using the following format:
        `<Column name>,<Threshold>:<Column name>,<Threshold>`

- `--no-default`
    Remove all default options. If set, PRSice will not set any defaults.

- `--pvalue` | `-p`

    Column header containing the p-value.
    The p-value information must be provided 

- `--se`

    Column header containing the standard error.
    Currently, this serves no purpose. 

- `--snp`

    Column header containing the SNP ID.
    This is required to allow SNP matching between the base and target file.

    !!! Note
        While it is possible to implement a feature to allow SNP matching purely based
        on the chromosome number and coordinate of a variant, the possibiliy of 
        flipping and multi-allelic input complicates the matter. Therefore this feature
        will not be implemented until an elegant solution can be provided.


- `--stat`

    Column header containing the summary statistic. If `--beta` is set, default
    as **BETA**. Otherwise, try and search for **OR** or **BETA** from the header of the base
    file. If both **OR** and **BETA** is presented in the header, PRSice will terminate. 
