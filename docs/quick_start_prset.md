# Background
A new feature of PRSice is the ability to perform gene set/pathway based analysis. This new feature is called PRSet.

!!! Important
    If you use the PRSet feature please cite our dedicated paper [PRSet: Pathway-based polygenic risk score analyses and software.](https://doi.org/10.1371/journal.pgen.1010624)

## Preparation
PRSet is based on [PRSice](quick_start.md), with additional input requirements

### Input
- **PRSice.R file**: A wrapper for the PRSice binary and for plotting
- **PRSice binary file**: Perform all analysis except plotting
- **Base data set**: GWAS summary results, which the PRS is based on
- **Target data set**: Raw genotype data of "target phenotype". Can be in the form of  [PLINK binary](https://www.cog-genomics.org/plink2/formats#bed) or [BGEN](http://www.well.ox.ac.uk/~gav/bgen_format/)

### PRSet Specific Input
- **Bed file(s)**: Bed file(s) containing region of genes within a gene set; or
- **MSigDB file**: File containing name of each gene sets and the ID of genes within the
gene set on each individual line. If MSigDB is provided, GTF file is required.
- **GTF file**: A file contain the genome boundary of each individual gene
- **SNP file**: A file containing SNPs constituting the gene set of interest. Can be in MSigDB (gmt) format or a file contain a single column of SNP IDs.

## Running PRSet

In most case, assuming the PRSice binary is located in `($HOME)/PRSice/bin/` and the working directory is `($HOME)/PRSice`, you can run PRSet with the following commands:

### With MSigDB data
Assuming a MSigDB file (*set.txt*) is [downloaded](http://software.broadinstitute.org/gsea/msigdb/) and a gene gtf file (gene.gtf) from [Ensemble](http://www.ensembl.org/index.html) is available, PRSet can then be performed using: 

``` bash hl_lines="7 8 9"
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

### With Bed Files
Alternatively, if a list of bed files are available, e.g. *A.bed,B.bed*, PRSet can be performed by running

``` bash hl_lines="7 8"
Rscript PRSice.R \
    --prsice ./bin/PRSice \
    --base TOY_BASE_GWAS.assoc \
    --target TOY_TARGET_DATA \
    --binary-target T \
    --thread 1 \
    --bed A.bed:SetA,B.bed \
    --multi-plot 10
```

!!! Note
    Both bed and GTF+MSigDB input can be used together

!!! Tips
    Name of the set will be the bed file name or can be provided using `--bed File:Name`

### With SNP Set
Finally, if you want to construct sets based on a list of SNPs, you can use `--snp-set`:

``` bash hl_lines="7 8"
Rscript PRSice.R \
    --prsice ./bin/PRSice \
    --base TOY_BASE_GWAS.assoc \
    --target TOY_TARGET_DATA \
    --binary-target T \
    --thread 1 \
    --snp-set A.snp:A,B.snp \
    --multi-plot 10
```

Two different format are allowed:

- SNP Set list format: A file containing a single column of SNP ID. Name of the set will be the file name or can be provided using `--snp-set File:Name`
- MSigDB format: Each row represent a single SNP set with the first column containing the name of the SNP set.
