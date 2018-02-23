
# Preparation
Before performing PRSice, you should perform quality control on your target samples. [See here](quick_start.md#quality-control-of-target-samples) for an example.

## Input
- **PRSice.R file:** A wrapper for the PRSice binary and for plotting
- **PRSice executable file:** Perform all analysis except plotting
- **Base data set:** GWAS summary results, which the PRS is based on
- **Target data set:** Raw genotype data of **target phenotype**.
Can be in the form of  [PLINK binary](https://www.cog-genomics.org/plink2/formats#bed) or [BGEN](http://www.well.ox.ac.uk/~gav/bgen_format/)

# Running PRSice
In most case, you can simply run PRSice using the following command, assuming your
PRSice executable is located in `($HOME)/PRSice/` and you are working in `($HOME)/PRSice`

!!! Note
    For window users, please use **Rscript.exe** instead of **Rscript**

## Binary Traits
For binary traits, you can use the following command (commands specific to binary traits are highlighted in yellow)
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
For quantitative traits, you can use  (commands specific to quantitative traits are highlighted in yellow)

``` bash hl_lines="6 7"
Rscript PRSice.R --dir . \
    --prsice ./PRSice \
    --base TOY_BASE_GWAS.assoc \
    --target TOY_TARGET_DATA \
    --thread 1 \
    --stat BETA \
    --binary-target F
```

!!! Note
    If you do not provide the PRSice command line with the type of Effect (`--stat`) or data type (`--binary-target`) then PRSice will try to work these out from the header of your base file:

    1. When *BETA* (case insensitive) is found in the header and `--stat` was not provided, `--beta` will be added to your command, and if `--binary-target` was not provided, `--binary-target F` will be added to your command 
   
    2. When *OR* (case insensitive) is found in the header and `--binary-target` was not provided, `--binary-target T` will be added to your command

    PRSice will detail all effective options in its log file where you can simply copy and paste it to get the same output

# Quality Control of Target Samples

You can perform quality control on the target samples using PLINK. 
A good starting point is (assume **_($target)_** is the prefix of your target binary file)

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

Then you can add `--keep ($target).qc.fam --extract ($target).qc.bim` to PRSice command to filter out
the samples and SNPs
