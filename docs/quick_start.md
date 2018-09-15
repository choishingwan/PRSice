
# Preparation
Before performing PRSice, quality control should be performed on the target samples. [See here](quick_start.md#quality-control-of-target-samples) for an example.

## Input
- **PRSice.R file:** A wrapper for the PRSice binary and for plotting
- **PRSice executable file:** Perform all analysis except plotting
- **Base data set:** GWAS summary results, which the PRS is based on
- **Target data set:** Raw genotype data of **target phenotype**.
Can be in the form of  [PLINK binary](https://www.cog-genomics.org/plink2/formats#bed) or [BGEN](http://www.well.ox.ac.uk/~gav/bgen_format/)

# Running PRSice
In most case, PRSice can simply be run using the following command, assuming the
PRSice executable is located in `($HOME)/PRSice/` and the working directory is `($HOME)/PRSice`

!!! Note
    For window users, please use **Rscript.exe** instead of **Rscript**

!!! Important
    Do not copy codes to Microsoft Word. Word has a tendency to change characters from codes into special characters that cannot be recognized by the terminal

## Binary Traits
For binary traits, the following command can be used (commands specific to binary traits are highlighted in yellow)
``` bash hl_lines="6 7"
Rscript PRSice.R --dir . \
    --prsice ./PRSice \
    --base TOY_BASE_GWAS.assoc \
    --target TOY_TARGET_DATA \
    --thread 1 \
    --stat OR \
    --binary-target T
```

## Quantitative Traits
For quantitative traits, the following can be used instead  (commands specific to quantitative traits are highlighted in yellow)

``` bash hl_lines="6 7 8"
Rscript PRSice.R --dir . \
    --prsice ./PRSice \
    --base TOY_BASE_GWAS.assoc \
    --target TOY_TARGET_DATA \
    --thread 1 \
    --stat BETA \
    --beta \
    --binary-target F
```

!!! Note
    If the type of Effect (`--stat`) or data type (`--binary-target`) were not specified, PRSice will try to determine these information based on the header of the base file:

    1. When *BETA* (case insensitive) is found in the header and `--stat` was not provided, `--beta` will be added to the command, and if `--binary-target` was not provided, `--binary-target F` will be added to the command 
   
    2. When *OR* (case insensitive) is found in the header and `--binary-target` was not provided, `--binary-target T` will be added to the command

    3. PRSice cannot determine if the type of effect / data type if the base file contains both *OR* and *BETA*

    PRSice will detail all effective options in its log file.

# Quality Control of Target Samples

Quality controls can be performed on the target samples using PLINK. 
A good starting point is (assume **_($target)_** is the prefix of the target binary file)

``` bash
plink --bfile ($target) \
    --maf 0.05 \
    --mind 0.1 \
    --geno 0.1 \
    --hwe 1e-6 \
    --make-just-bim \
    --make-just-fam \
    --out ($target).qc
```

Then, `--keep ($target).qc.fam --extract ($target).qc.bim` can be added to the PRSice command to filter out
the samples and SNPs
