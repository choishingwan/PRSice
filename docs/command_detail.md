This page contains all command available in PRSice. 

!!! tips
    
    When constructing new parameters, we follow the following rule: if the command has effect on any file that is not the target, 
    it will have a prefix of the file name. For example, `--base-info` applies INFO score filtering on the base file, `--ld-info` perform INFO score filtering on the LD reference file and `--info` applies the INFO score filtering on the target file. 


# Base File
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
    of BETA. If set, PRSice assume the statistic is in the form BETA. 
    Mutually exclusive from `--or`

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

- `--base-info`

    Base INFO score filtering. Format should be `<Column name>,<Threshold>`.
    SNPs with info score less than `<Threshold>` will be ignored.
    It is useful to perform INFO score filtering to remove SNPs
    with low imputation confidence score. By default, PRSice will search for
    the **INFO** column in your base file and perform info score filtering with
    threshold of **0.9**. You can disable this behaviour by using `--no-default`
    
- `--base-maf`

    Base minor allele frequency (MAF) filtering.
    Format should be `<Column name>,<Threshold>`.
    SNPs with MAF less than `<Threshold>` will be ignored.
    Additional column can be provided (e.g. different filtering threshold
    for case and control), using the following format:
        `<Column name>,<Threshold>:<Column name>,<Threshold>`


- `--no-default`
    Remove all default options. If set, PRSice will not set any defaults.

- `--or`

    This flag is used to indicate if the test statistic is in the form
    of odd ratios. If set, PRSice assume the statistic is in the form OR. 
    Mutually exclusive from `--beta`

- `--pvalue` | `-p`

    Column header containing the p-value.
    The p-value information must be provided 

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
    to **BETA**; likewise, if `--or` is set, default to **OR**. 
    Otherwise, will try and search for **OR** or **BETA** from the header of the base
    file. If both **OR** and **BETA** is presented in the header, PRSice will terminate. 

# Target File

- `--binary-target`

    Indicate whether the target phenotype is binary or not.
    Either **T** or **F** should be provided where **T** represent a binary phenotype.
    For multiple phenotypes, the input should be separated by comma without space.

    Default: **F** if `--beta` is set and **T** if `--or` is set


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

    Filter SNPs based on minor allele frequency (MAF). MAF is calculated using only the founder samples

    !!! Note

        When perform MAF filtering on dosage data, the MAF 
        is calculated using the hard-coded genotype


- `--nonfounders`

    By default, PRSice will exclude all non-founders from the analysis.
    When this flag is set, non-founders will be included in the
    regression model but will still be excluded from LD estimation.

- `--pheno` | `-f`

    Tab or space delimited phenotype file containing the phenotype(s).
    First column must be FID of the samples and
    the second column must be IID of the samples.
    When `--ignore-fid` is set, first column must
    be the IID of the samples. Must contain a
    header if `--pheno-col` is specified

- `--pheno-col`

    Headers of phenotypes to be included from the phenotype file.
    When multiple phenotypes are provided, the phenotype name will be
    used as part of the file output prefix

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
    `--target <prefix>,<fam or sample file>`

- `--target-list`

    File containing prefix of target genotype files. 
    Similar to `--target` but allow for more flexibility.
    A separate fam/sample file can be specified by
    `--target-list <list-file>,<fam or sample file>`

- `--type`

    File type of the target file. Support bed (binary plink) and bgen format. Default: bed

# Dosage
- `--allow-inter`

    Allow the generate of intermediate file. This will
    speed up PRSice when using dosage data as clumping
    reference and for hard coding PRS calculation

- `--dose-thres`

    Translate any SNPs with highest genotype probability less than this threshold to missing call. 
    For example, with `--dose-thres 0.9`, sample with genotype probability of $P(0/0)=0.2$, $P(0/1)=0.52$, $P(1/1)=0.28$ will be set to missing

- `--hard-thres`

    A hardcall is saved when the distance to the nearest 
    hardcall is less than the hardcall threshold.
    Otherwise a missing code is saved. Default is: 0.1

    The distance ($D$) to the nearest hardcall is calculated as:

    $$
        P(Ref) = 2 \times P(HomRef) + P(Het) \\
        P(Alt) = 2 \times P(HomAlt) + P(Het) \\
        D = 0.5 \times \left(|P\left(Ref\right)- round\left(P\left(Ref\right)\right)| + |P\left(Alt\right)- round\left(P\left(Alt\right)\right)|\right)
    $$

    !!! Note

        If dosage data is used as a LD reference, it will always be
        hard coded to calculate the LD

        Default: 0.9

- `--hard`

    When set, will use hard thresholding instead of dosage for PRS construction.
    Default is to use dosage.

# Clumping
- `--clump-kb`
    The distance for clumping in kb.
    For example, if `--clump-kb 250` is provided, PRSice will clump any SNPs that is 
    within 250kb to **both** end of the index SNP (therefore a 500kb window with the index SNP at the center).
    Now also support distance with a unit. e.g. `--clump-kb 1M` is a valid input.
    Default: 1M

- `--clump-r2`

    The r^2^ threshold for clumping. Default: 0.1

- `--clump-p`

    The p-value threshold use for clumping. Default: 1.

- `--ld` | `-L`

    LD reference file. Use for estimation of LD during clumping.
    If not provided, will use the post-filtered target genotype
    for LD calculation. Support multiple chromosome input.
    Please see [`--target`](command_detail.md#target-file) for more information.
    When the target sample is small (e.g. < 500) and
    external panel of the same population is available (e.g. 1000 genome),
    an external reference panel might be used to improve the LD estimation for clumping.

- `--ld-dose-thres`

    Translate any SNPs with highest genotype probability less than this threshold to missing call. 
    For example, with `--ld-dose-thres 0.9`, sample with genotype probability of $P(0/0)=0.2$, $P(0/1)=0.52$, $P(1/1)=0.28$ will be set to missing

- `--ld-geno`

    Filter SNPs based on genotype missingness. Must be a value
    between *0.0* and *1.0*.

- `--ld-info`

    Filter SNPs based on info score. Only used for imputed LD reference.
    The INFO score is calculated as the MaCH imputation r-squared value, 
    represented by the following pseudo code

```
    m=Mean of expected genotype
    v=variance of expected genotype
    p=m/2
    p_a = 2p(1-p)
    INFO = v/p_a
```

- `--ld-hard-thres`

    A hardcall is saved when the distance to the nearest 
    hardcall is less than the hardcall threshold.
    Otherwise a missing code is saved. Default is: 0.1

    The distance ($D$) to the nearest hardcall is calculated as:

    $$
        P(Ref) = 2 \times P(HomRef) + P(Het) \\
        P(Alt) = 2 \times P(HomAlt) + P(Het) \\
        D = 0.5 \times \left(|P\left(Ref\right)- round\left(P\left(Ref\right)\right)| + |P\left(Alt\right)- round\left(P\left(Alt\right)\right)|\right)
    $$

    !!! Note

        If dosage data is used as a LD reference, it will always be
        hard coded to calculate the LD

        Default: 0.9

- `--ld-keep`

    File containing the sample(s) to be extracted from
    the LD reference file. First column should be FID and
    the second column should be IID. If `--ignore-fid` is
    set, first column should be IID.
    Mutually exclusive from `--ld-remove`. No effect if `--ld` was not provided

- `--ld-list`

    File containing prefix of multiple LD reference files. 
    Similar to --ld but allow more flexibility. 
    A separate fam/sample file can be specified by
    `--ld-list <list-file>,<fam or sample file>`

- `--ld-maf`

    Filter SNPs based on minor allele frequency (MAF)

    !!! Note

        When perform MAF filtering on dosage data, MAF 
        is calculated using the hard-coded genotype

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

    When set, PRSice will not perform clumping. This is useful
    a pre-clumped list of SNPs is available.

- `--pearson`

    Use Pearson Correlation instead of maximum likelihood haplotype frequency estimates for LD calculation

- `--proxy`

    Proxy threshold for index SNP to be considered
    as part of the region represented by the clumped
    SNP(s). e.g. `--proxy 0.8` means the index SNP will
    represent region of any clumped SNP(s) that has
    r^2^=0.8 even if the index SNP does not physically
    locate within the region


# Covariate
- `--cov` | `-C`
    Covariate file. First column should be FID and
    the second column should be IID. If `--ignore-fid`
    is set, first column should be IID

- `--cov-col` | `-c`
    Header of covariates. If not provided, will use
    all variables in the covariate file. By adding
    `@` in front of the string, any numbers within `[`
    and `]` will be parsed. E.g. `@PC[1-3]` will be
    read as **PC1,PC2,PC3**. Discontinuous input are also
    supported: `@cov[1.3-5]` will be parsed as
    **cov1,cov3,cov4,cov5**
    
- `--cov-factor`

    Header of categorical covariate(s). Dummy variable 
    will be automatically generated. Any items in
    `--cov-factor` must also be found in `--cov-col`
    Also accept continuous input (start with `@`).

#P-value Thresholding
- `--bar-levels`

    Level of barchart to be plotted. When `--fastscore`
    is set, PRSice will only calculate the PRS for
    threshold within the bar level. Levels should be
    comma separated without space

- `--fastscore`

    Only calculate threshold stated in `--bar-levels`

- `--no-full`

    By default, PRSice will include the full model,
    i.e. p-value threshold = 1. Setting this flag will
    disable that behaviour

- `--interval` | `-i`

    The step size of the threshold. Default: 5e-05

- `--lower` | `-l`

    The starting p-value threshold. Default: 5e-08

- `--model`

    Genetic model use for regression. The genetic
    encoding is based on the base data where the
    encoding represent number of the effective allele
    Available models include:

    - `add` - Additive model, code as 0/1/2 (default)
    - `dom` - Dominant model, code as 0/1/1
    - `rec` - Recessive model, code as 0/0/1
    - `het` - Heterozygous only model, code as 0/1/0

- `--missing`

    Method to handle missing genotypes. 
    Available methods include:
    - `MEAN_IMPUTE` - Missing genotypes contribute an amount
    proportional to imputed allele frequency (default)
    - `SET_ZERO` - To throw out missing observations instead
    - `CENTER` - shift all scores to mean zero.

- `--no-regress`            

    Do not perform the regression analysis and simply
    output all PRS.

- `--score`

    Method to calculate the polygenic score.
    Available methods include:

    - `avg` - Take the average effect size (default)
    - `std` - Standardize the effect size
    - `sum` - Direct summation of the effect size

- `--upper` | `-u`

    The final p-value threshold. Default: 0.5

# PRSet
- `--background`

    
    String to indicate a background file. This string
    should have the format of Name:Type where type can be

    - bed   - 0-based range with 3 column. Chr Start End
    - range - 1-based range with 3 column. Chr Start End
    - gene  - A file contain a column of gene name

    As the name suggest, the background file inform PRSet of 
    the background signal to be used for competitive p-value calculation.
    
    When a background file is not provided, PRSet will construct the background
    using the GTF file. However, if both the background and the GTF file isn't provided, 
    PRSet cannot perform the set base permutation. In this case, you can use `--full-back` to 
    indicate that you'd like to use the whole genome as the background set
  
- `--bed` | `-B`

    Bed file containing the selected regions.
    Name of bed file will be used as the region
    identifier. 
    
    !!! warning

        Bed file is 0-based. 

- `--feature`

    Feature(s) to be included from the gtf file.

    Default: *exon,CDS,gene,protein_coding*

- `--full-back`

    Use the whole genome as background set for competitive p-value calculation

- `--gtf` | `-g`

    GTF file containing gene boundaries. Required when `--msigdb` is used

    !!! tip

        Human Genome build GRCh38 can be downloaded from [here](ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens).

- `--msigdb` | `-m`

    MSigDB file containing the pathway information. Require the gtf file.
    The GMT file format used by MSigDB is a simple tab/space delimited text file
    where each line correspond to a single gene set following by Gene IDs:
    ```
        [Set A] [Gene 1] [Gene 2] ...
        [Set B] [Gene 1] [Gene 2] ...
    ```

    !!! tip

        Curated MSigDB files can be downloaded from [here](http://software.broadinstitute.org/gsea/msigdb/) after registration in [here](http://software.broadinstitute.org/gsea/login.jsp;jsessionid=EEFB5FCE8B9B285B2F789B46B388A647#msigdb)

- `--set-perm`               

    The number of set base permutation to perform. 
    This is only used for calculating the competitive p-value. 
    10,000 permutation nshould generally be enough. 

- `--snp-set`               

    Provide gene sets using SNP ID. Two different format is allowed:

    - SNP Set list format: A file containing a single column of SNP ID. Name of the set will
    be the file name or can be provided using `--snp-set File:Name`
    - MSigDB format: Each row represent a single SNP set with the first column containing the name of the SNP set. 
    
- `--wind-3`

    Add N base(s) to the 3' region of each gene regions. Unit suffix are allowed e.g. `--wind-3 1M`
    
- `--wind-5`

    Add N base(s) to the 5' region of each gene regions. Unit suffix are allowed e.g. `--wind-5 1M`

# R specific commands
- `--prsice`

    Location of the PRSice executable.

- `--dir`

    Location to install require R packages. Only require if the required packages are not installed.
    We require the following packages: `optparse`, `method`, `tools`, `ggplot2`, `data.table`, `grDevices`, `RColorBrewer`
    
# Plotting 

- `--bar-col-high`

   Colour of the most predicting threshold. Default: `firebrick`

- `--bar-col-low`

  Colour of the poorest predicting threshold. Default: `dodgerblue`

- `--bar-col-p`

  When set, will change the colour of bar to p-value threshold instead of
  the p-value from the association with phenotype

- `--bar-palatte`

  Colour palatte to be used for bar plotting when `--bar_col_p` is set. Default: `YlOrRd`

- `--device`

    Select different plotting devices. You can choose
    any plotting devices supported by base R. Default: png

- `--multi-plot`

  Plot the top N target phenotypes / gene sets in a summary plot


- `--plot`

  When set, will only perform plotting using existing PRSice result files. 
  All other parameters are still required such that PRSice
  can correctly locate the required input files for plotting.

- `--plot-set`

  The default behaviour of PRSet is to plot the bar-chart, high-resolution plot and
  quantile plot of the "Base" gene set, which consider
  all SNPs within the genome. By using the `--plot-set` option, you can plot the
  specific set of interest.

- `--quantile` | `-q`

    Number of quantiles to plot.
    No quantile plot will be generated when this is not provided.

-  `--quant-break`

    Parameter to indicate an uneven distribution of quantile. Values represent
    the upperbound of each quantile group. 

    e.g. With `--quantile 10 --quant-break 1,5,10`, the quantiles will be grouped into

    > $0\lt Q \le 1$, $1\lt Q \le 5$, $5\lt Q \le 10$

    !!! Note
    
        To use `--quant-break`, you must set the correct amount of quantiles. For example, 
        if the largest value in `--quant-break` is 100, then you must use `--quantile 100`

- `--quant-extract` | `-e`

  File containing sample ID to be plot on a separated
  quantile e.g. extra quantile containing only
  schizophrenia samples. Must contain IID. Should
  contain FID if `--ignore-fid` isn't set.

!!! note

    This will only work if the base and target has a different
    phenotype or if the target phenotype is quantitative


- `--quant-ref`

  Reference quantile for quantile plot. Default is number of quantiles divided by 2

  Or in the event where `--quant-break` is used, represent the upper bound of the 
  reference quantile

- `--scatter-r2`

  When set, will change the y-axis of the high resolution scatter plot to R2 instead


# Miscellaneous
- `--all-score`

    Output PRS for ALL threshold.

    !!! warning
        This will generate a huge file

- `--enable-mmap`
           
    Enable memory mapping. This will provide a small speed boost if large amount of memory
    is available and if all genotypes were stored in the same file

    !!! Warning
        If you don't have enough memory, and if your genotypes were stored in different files, 
        memory mapping will actually slow down PRSice. 

- `--exclude`

    File contains SNPs to be excluded from the analysis.
    Mutually exclusive from `--extract`

- `--extract`

    File contains SNPs to be included in the analysis.
    Mutually exclusive from `--exclude`

- `--id-delim`

    Delimiter used to concatinate FID and IID when performing ID matching. 
    Especially useful for BGEN file processing

- `--ignore-fid`

    Ignore FID for all input. When this is set,
    first column of all file will be assume to
    be IID instead of FID

    
- `--keep-ambig`

    Keep ambiguous SNPs. Only use this option
    whe base and target has the same A1 and A2 alleles

- `--logit-perm`

    When performing permutation on binary phenotypes, 
    use logistic regression instead of linear regression. 
    This will substantially slow down PRSice.

    !!! note

        One problem with using `--logit-perm` is that some of
        the permuted phenotype might be suffer from perfect
        separation. This leads to the GLM logistic model not
        being able to be converge (thus terminating PRSice).

        If you encounter such problem, you might want to exclude
        the `--logit-perm` option. In most case, the p-value of the
        linear model should be similar to the logistic model

- `--memory`
    
    Maximum memory usage allowed. PRSice will try its best to honor this setting. 
    For example, `--memory 10Gb` will restrict PRSice to use no more than 10Gb of memory.  
    However, as we are not using memory pool like PLINK, it is possible for PRSice 
    to use more than the allowed amount. PRSice will mainly check the memory usage when:

    - Perform Clumping
    - Perform permutation analysis
    - Perform set-based permutation
 
- `--non-cumulate`
    
    Calculate non-cumulative PRS. PRS will be reset
    to 0 for each new P-value threshold instead ofadding up

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
    use value larger than or equal to 10,000

    !!! note

        When permutation is required, PRSice will perform the following
        operation

        1. Perform normal PRSice across all thresholds and obtain p-value of the
        most significant threshold
        2. Repeat PRSice analysis *N* times with permuted phenotype. Count the
        number of time where the p-value of the most significant threshold for
        the permuted

- `--print-snp`

    Print all SNPs that remains in the analysis after clumping is performed. For PRSet, `Y` indicate the SNPs
    falls within the gene set of interest and `N` otherwise. If only PRSice is performed, a single "gene set" called
    "Base" will be indicated with all entries marked as `Y`

- `--seed` | `-s`

    Seed used for permutation. If not provided,
    system time will be used as seed. This will
    allow the same results to be generated when
    the same seed and input is used

- `--thread` | `-n`

    Number of thread use

    !!! tip

        Maximum number of thread can be specified by using `--thread max`

    !!! note

        PRSice will limit the maximum number of thread used to the number of core available on the system as detected by PRSice.

- `--x-range`               
    Range of SNPs to be excluded from the whole
    analysis. It can either be a single bed file
    or a comma seperated list of range. Range must
    be in the format of *chr:start-end* or *chr:coordinate*

- `--help` | `-h`

    Display the help messages
