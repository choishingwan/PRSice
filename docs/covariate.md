# Introduction
Systematic ancestry difference among samples can cause spurious associations and
can lead to invalid PRS results. Similarly, power of the PRS analysis can be 
improved by controlling for other confounders. To account for that,
PRSice allow the incorporation of covariates into the analysis.

!!! note

    When large number of covariates are included in the model, missing
    data might poses a problem. PRSice will automatically exclude any
    samples with missing covariates from the regression model but will
    still calculate the PRS for that sample. The **In_Regression** column
    in the best score output is used to indicate whether the sample is
    included in the regression model (*Yes* for included; *No* for excluded)

# Comamnds
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

- `--cov-file` | `-C`
    Covariate file. First column should be FID and
    the second column should be IID. If `--ignore-fid`
    is set, first column should be IID
