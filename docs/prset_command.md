# Introduction

!!! important

    It is vital that the human genome build is the same for the GTF file,
    bed files, target file and the base file.
    Otherwise the coordinates of the SNPs can be wrong and PRSice will
    not be able to correctly assign the gene membership, leading to
    invalid results.

# Command
- `--bed` | `-B`

    Bed file containing the selected regions.
    Name of bed file will be used as the region
    identifier. WARNING: Bed file is 0-based

- `--feature`

    Feature(s) to be included from the gtf file.

    Default: *exon,CDS,gene,protein_coding*

- `--gtf` | `-g`

    GTF file containing gene boundaries. Required when `--msigdb` is used

    !!! tip

        Human Genome build GRCh38 can be downloaded from [here](ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens).

- `--msigdb` | `-m`

    MSIGDB file containing the pathway information. Require the gtf file.
    The MSigDB file format is a simple tab/space delimited text file
    where each line correspond to a single gene set following by Gene IDs:
    ```
        [Set A] [Gene 1] [Gene 2] ...
        [Set B] [Gene 1] [Gene 2] ...
    ```

    !!! tip

        Curated MSigDB files can be downloaded from [here](http://software.broadinstitute.org/gsea/msigdb/) after registration in [here](http://software.broadinstitute.org/gsea/login.jsp;jsessionid=EEFB5FCE8B9B285B2F789B46B388A647#msigdb)
