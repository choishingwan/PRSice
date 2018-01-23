# Introduction
Base (i.e. GWAS) data must be provided as a whitespace delimited file containing association analysis results for SNPs on the base phenotype.
PRSice has no problem reading in a gzipped base file.

If you are using output from PLINK, then please make sure there is a column for the effect allele (A1) and specify it with `--A1` option.
On the other hand, if you are using another form of output, then you can provide the column headers using the `--chr`, `--A1`, `--A2`, `--stat`, `--snp`, `--bp`, `--se`, `--pvalue` options

!!! important

    For PRSice to run, the base file must contain the effect allele (`--A1`), effect size estimates (`--stat`), p-value for association (`--pvalue`), and the SNP id (`--snp`).

If the input file does not contain a column header, then you can specify the column using their index (start counting from 0) with the `--index` flag.

e.g. if the first column of the file contains the p-value information, you can use:
``
--pvalue 0 --index
``

Below is an example input file:

|SNP|CHR|BP|A1|A2|OR|SE|P|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|rs3094315|1|752566|A|G|0.9912|0.0229|0.7009|
|rs3131972|1|752721|A|G|1.007|0.0228|0.769|
|rs3131971|1|752894|T|C|1.003|0.0232|0.8962|

and the parameters can either be

``
--snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat OR --se SE --pvalue P
``

or

``
--snp 0 --chr 1 --bp 2 --A1 3 --A2 4 --stat 5 --se 6 --pvalue 7 --index
``

Strand flips are automatically detected and accounted for

If an imputation info score or the minor allele frequencies (MAF) are also included in the file, you can use
`--info-base <Info Name>,<Info Threshold>` and `--maf-base <MAF Name>,<MAF Threshold>` to filter
SNPs based on their INFO score and MAF respectively.

For binary trait base file, you can also filter the MAF according the the case and control separately using

```
--maf-base <MAF Name for Control>,<MAF Threshold for Control>:<MAF Name for Case>,<MAF Threshold for Case>
```

By default, PRSice will look for the following column names automatically from your base file header if `--index` was not provided:
> CHR, BP, A1, A2, SNP, P, INFO, SE (case sensitive) and OR / BETA (case insensitive)

If you don't want the specific column to be included e.g. You don't want to perform info score filtering but your file contains an
INFO column, you can use `--no-default` to disable all the defaults of PRSice.
PRSice will ignore any columns that were not found in your file (e.g. If you specified `--A2 B` but none of your column header is *B*,  then PRSice will treat it as if no *A2* information is presented)



# Commands
- `--A1`

    Column header containing the **effective allele**.
    As there is not a standardized name, please make sure you have
    provided the correct column. Otherwise, the effect will be flipped.

- `--A2`

    Column header containing **non-effective allele**.
    As there is not a standardized name, please make sure you have
    provided the correct column. Otherwise, the effect will be flipped.

- `--base` | `-b`

    Base (i.e. GWAS) association file. This is a whitespace delimited file
    containing association results for SNPs on the base phenotype.
    This file can be *gzipped*.
    For PRSice to run, the base file must contain the effect allele
    (`--A1`), effect size estimates (`--stat`), p-value for association
    (`--pvalue`), and the SNP ID (`--snp`).

- `--beta`

    This flag is used to indicate if the test statistic is in the form
    of BETA. If set, test statistic is assume to be in the form of BETA

- `--bp`

    Column header containing the coordinate of SNPs. When this is provided,
    the coordinate of the SNPs will be scrutinized between the base and target file.
    SNPs with mismatched coordinate will be automatically excluded.

- `--chr`

    Column header containing the chromosome information. When this is provided,
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
    An additional column can also be added (e.g. different filtering threshold
    for case and control), using the following format:
        `<Column name>,<Threshold>:<Column name>,<Threshold>`

- `-no-default`
    Remove all default options. If set, PRSice will not set any defaults
    and you will have to ensure all required columns are provided.
    (`--snp`, `--stat`, `--A1`, `--pvalue`)

- `--pvalue` | `-p`

    Column header containing the p-value.
    The p-value information is required for clumping and
    this information must be provided under all circumstances.

- `--se`

    Column header containing the standard error.
    Currently we do not use the standard error in any of the analysis.
    Thus this column is not required.

- `--snp`

    Column header containing the SNP ID.
    This is required such that we can match SNPs between the base and target file.


- `--stat`

    Column header containing the summary statistic. If `--beta` is set, default
    as **BETA**. Otherwise, try and search for **OR** or **BETA** from the header of the base
    file. The summary statistic is required for the construction of PRS.
