# Introduction
A PRS for an individual is a summation of their genotypes at variants
genome-wide, weighted by effect sizes on a trait of interest.
Effect sizes are typically estimated from published GWAS results,
and only variants exceeding a P-value threshold, PT, are included.
Since even large GWAS achieve only marginal evidence for association
for many causal variants, PRS are usually calculated at a set of
P-value thresholds, e.g.  PT=1×10^−5^,1×10^−4^,…,0.05,0.1,…,0.5

PRSice will automatically calculate the PRS for different p-value thresholds
and perform a regression to test the level of association of the PRS with
the target phenotype. This allow the identification of PRS that "best" predicts
the phenotype and can be used for downstream analysis.

!!! Note

    The p-value thresholds are inclusive. That is, For a p-value threshold of
    0.5, all SNPs with p-value of 0.5 will also be included in this threshold.

# Command
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

    The step size of the threshold. Default: 0.00005

- `--lower` | `-l`

    The starting p-value threshold. Default: 0.0001

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
    - `no_mean_imputation` - To throw out missing observations instead
    - `center` - shift all scores to mean zero.

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
