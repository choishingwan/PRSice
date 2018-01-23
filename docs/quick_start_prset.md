# Background
A new feature of PRSice is the ability to perform set base/pathway based analysis. This new feature is called PRSet.

Paper currently under-preparation


!!! Important
    This feature is currently under active development. Most notably, the current
    statistical testing is only a simple regression. A better statistic method
    will soon be implemented. Meanwhile, we recommend using [MAGMA](https://ctg.cncr.nl/software/magma)
    to identify significant gene sets and use these gene sets as an input to PRSet.

# Preparation
PRSet is mainly based on [PRSice](quick_start.md). To perform PRSet,
you will need the following inputs

## Input
- **PRSice.R file**: A wrapper for the PRSice binary and for plotting
- **PRSice binary file**: Perform all analysis except plotting
- **Base data set**: GWAS summary results, which the PRS is based on
- **Target data set**: Raw genotype data of "target phenotype". Can be in the form of  [PLINK binary](https://www.cog-genomics.org/plink2/formats#bed) or [BGEN](http://www.well.ox.ac.uk/~gav/bgen_format/)
## PRSet Specific Input
- **Bed file(s)**: Bed file(s) containing region of genes within a gene set; or
- **MSigDB file**: File containing name of each gene sets and the ID of genes within the
gene set on each individual line. If MSigDB is provided, GTF file is required.
- **GTF file**: A file contain the genome boundary of each individual gene

!!! Note

    Currently, it is recommended to perform MAGMA to identify significant gene sets and use those as an input of PRSet

# Running PRSet

In most case, you can simply run PRSice using the following command, assuming your
PRSice binary is located in `($HOME)/PRSice/bin/` and you are working in `($HOME)/PRSice`

## With MSigDB data
Assuming you have [downloaded](http://software.broadinstitute.org/gsea/msigdb/) a MSigDB file (*set.txt*) and a gene gtf file (gene.gtf) from [Ensemble](http://www.ensembl.org/index.html), then you can run PRSet with the following command

```
Rscript PRSice.R \
    --prsice ./bin/PRSice \
    --base TOY_BASE_GWAS.assoc \
    --target TOY_TARGET_DATA \
    --binary-target T \
    --thread 1 \
    --gtf gene.gtf \
    --msigdb set.txt \
    --multi-plot 10
```

This will perform PRSet analysis and generate the multi-set plot with the top 10 gene sets

## With Bed Files
Alternatively, if you have a list of bed files e.g. *A.bed,B.bed* as an input, you can run PRSet as

```
Rscript PRSice.R \
    --prsice ./bin/PRSice \
    --base TOY_BASE_GWAS.assoc \
    --target TOY_TARGET_DATA \
    --binary-target T \
    --thread 1 \
    --bed A.bed,B.bed \
    --multi-plot 10
```

!!! Note
    You can use both bed and GTF+MSigDB input together*
