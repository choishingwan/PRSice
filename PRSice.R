#!/usr/bin/env Rscript
# Here is the guide to this protentially long R code
# To go to each section, just search for the corresponding header as stated here
# The code structure are as follow
# The easiest way will be to use RStudio and go to the corresponding section by
# selecting it on the bottom left corner of the script console
#
# INSTALL_PACKAGE
# - Contains functions responsible for installing all required packages
#
# COMMAND_FUNC
# - Functions required for command line argument parsing
#
# COMMAD_BUILD
# - Building the command line parser using argparser
#
# CALL_PRSICE
# - call the cpp prsice
#
# PLOTTING
# - Here contains all the function for plotting
# - quantile_plot: plotting the quantile plots
# - run_plot: The function used for calling different plotting functions
#
# CALL PLOTTING FUNCTION
# - Process the input names and call the actual plotting function
#
# Environment stuff: This will allow us to locate the cpp file correctly
#
# For window users, the reason why PRSice doesn't work with window is due to
# area flagged with WINDOW PEOPLE. If you can go around that, then you can
# make PRSice work (though I have not debug PRSice executable in window)


# Remove annoying messages ------------------------------------------------
#options(error = quote({dump.frames(to.file=TRUE); q()}))
In_Regression <-
    DEC <-
    Coef <-
    CI.L <-
    CI.U <-
    Group <-
    Threshold <-
    R2 <-
    print.p <- R <- P <- value <- Phenotype <- Set <- PRS.R2 <- LCI <- UCI <- quant.ref <- NULL

# Help Messages --------------------------------------
help_message <-
"usage: Rscript PRSice.R [options] <-b base_file> <-t target_file> <--prsice prsice_location>\n
\nRequired:\n
    --prsice                Location of the PRSice binary\n
    --dir                   Location to install ggplot. Only require if ggplot\n
                            is not installed\n
\nBase File:\n
    --A1                    Column header containing allele 1 (effective allele)\n
                            Default: A1\n
    --A2                    Column header containing allele 2 (non-effective allele)\n
                            Default: A2\n
    --base          | -b    Base association file\n
    --beta                  Whether the test statistic is in the form of \n
                            BETA or OR. If set, test statistic is assume\n
                            to be in the form of BETA.\n
    --bp                    Column header containing the SNP coordinate\n
                            Default: BP\n
    --chr                   Column header containing the chromosome\n
                            Default: CHR\n
    --index                 If set, assume the INDEX instead of NAME  for\n
                            the corresponding columns are provided. Index\n
                            should be 0-based (start counting from 0)\n
    --info-base             Base INFO score filtering. Format should be\n
                            <Column name>,<Threshold>. SNPs with info \n
                            score less than <Threshold> will be ignored\n
                            Column name default: INFO\n
                            Threshold default: 0.9\n
    --maf-base              Base MAF filtering. Format should be\n
                            <Column name>,<Threshold>. SNPs with maf\n
                            less than <Threshold> will be ignored. An\n
                            additional column can also be added (e.g.\n
                            also filter MAF for cases), using the\n
                            following format:\n
                            <Column name>,<Threshold>:<Column name>,<Threshold>\n
    --no-default            Remove all default options. If set, PRSice\n
                            will not set any default column name and you\n
                            will have to ensure all required columns are\n
                            provided. (--snp, --stat, --A1, --pvalue)\n
    --pvalue        | -p    Column header containing the p-value\n
                            Default: P\n
    --se                    Column header containing the standard error\n
                            Default: SE\n
    --snp                   Column header containing the SNP ID\n
                            Default: SNP\n
    --stat                  Column header containing the summary statistic\n
                            If --beta is set, default as BETA. Otherwise,\n
                            will search for OR or BETA from the header\n
                            of the base file\n
\nTarget File:\n
    --binary-target         Indicate whether the target phenotype\n
                            is binary or not. Either T or F should be\n
                            provided where T represent a binary phenotype.\n
                            For multiple phenotypes, the input should be\n
                            separated by comma without space. \n
                            Default: T if --beta and F if --beta is not\n
    --geno                  Filter SNPs based on gentype missingness\n
    --info                  Filter SNPs based on info score. Only used\n
                            for imputed target\n
    --keep                  File containing the sample(s) to be extracted from\n
                            the target file. First column should be FID and\n
                            the second column should be IID. If --ignore-fid is\n
                            set, first column should be IID\n
                            Mutually exclusive from --remove\n
    --maf                   Filter SNPs based on minor allele frequency (MAF)\n
    --nonfounders           Keep the nonfounders in the analysis\n
                            Note: They will still be excluded from LD calculation\n
    --pheno-col     | -F    Headers of phenotypes to be included from the\n
                            phenotype file\n
    --pheno-file    | -f    Phenotype file containing the phenotype(s).\n
                            First column must be FID of the samples and\n
                            the second column must be IID of the samples.\n
                            When --ignore-fid is set, first column must\n
                            be the IID of the samples.\n
                            Must contain a header if --pheno-col is\n
                            specified\n
    --prevalence    | -k    Prevalence of all binary trait. If provided\n
                            will adjust the ascertainment bias of the R2.\n
                            Note that when multiple binary trait is found,\n
                            prevalence information must be provided for\n
                            all of them (Either adjust all binary traits,\n
                            or don't adjust at all)\n
    --remove                File containing the sample(s) to be removed from\n
                            the target file. First column should be FID and\n
                            the second column should be IID. If --ignore-fid is\n
                            set, first column should be IID\n
                            Mutually exclusive from --keep\n
    --target        | -t    Target genotype file. Currently support\n
                            both BGEN and binary PLINK format. For \n
                            multiple chromosome input, simply substitute\n
                            the chromosome number with #. PRSice will\n
                            automatically replace # with 1-22\n
                            For binary plink format, you can also specify\n
                            a seperate fam file by <prefix>,<fam file>\n
    --target-list           File containing prefix of target genotype\n
                            files. Similar to --target but allow more \n
                            flexibility. Do not support external fam file\n
                            at the moment\n
    --type                  File type of the target file. Support bed \n
                            (binary plink) and bgen format. Default: bed\n
\nDosage:\n
    --allow-inter           Allow the generate of intermediate file. This will\n
                            speed up PRSice when using dosage data as clumping\n
                            reference and for hard coding PRS calculation\n
    --hard-thres            Hard threshold for dosage data. Any call less than\n
                            this will be treated as missing. Note that if dosage\n
                            data is used as a LD reference, it will always be\n
                            hard coded to calculate the LD\n
    --hard                  Use hard coding instead of dosage for PRS construction.\n
                            Default is to use dosage instead of hard coding\n
\nClumping:\n
    --clump-kb              The distance for clumping in kb\n
                            Default: 250 \n
    --clump-r2              The R2 threshold for clumping\n
                            Default: 0.1 \n
    --clump-p               The p-value threshold use for clumping.\n
                            Default: 1.0 \n
    --ld            | -L    LD reference file. Use for LD calculation. If not\n
                            provided, will use the post-filtered target genotype\n
                            for LD calculation. Support multiple chromosome input\n
                            Please see --target for more information\n
    --ld-list               File containing prefix of LD reference files.\n
                            Similar to --ld but allow more \n
                            flexibility. Do not support external fam file\n
                            at the moment\n
    --ld-geno               Filter SNPs based on genotype missingness\n
    --ld-info               Filter SNPs based on info score. Only used\n
                            for imputed LD reference\n
    --ld-hard-thres         Hard threshold for dosage data. Any call less than\n
                            this will be treated as missing.\n
                            Default: 0.9\n
    --ld-keep               File containing the sample(s) to be extracted from\n
                            the LD reference file. First column should be FID and\n
                            the second column should be IID. If --ignore-fid is\n
                            set, first column should be IID\n
                            Mutually exclusive from --ld-remove\n
                            No effect if --ld was not provided\n
    --ld-maf                Filter SNPs based on minor allele frequency\n
    --ld-remove             File containing the sample(s) to be removed from\n
                            the LD reference file. First column should be FID and\n
                            the second column should be IID. If --ignore-fid is\n
                            set, first column should be IID\n
                            Mutually exclusive from --ld-keep\n
    --ld-type               File type of the LD file. Support bed (binary plink)\n
                            and bgen format. Default: bed\n
    --no-clump              Stop PRSice from performing clumping\n
    --proxy                 Proxy threshold for index SNP to be considered\n
                            as part of the region represented by the clumped\n
                            SNP(s). e.g. --proxy 0.8 means the index SNP will\n
                            represent region of any clumped SNP(s) that has a\n
                            R2>=0.8 even if the index SNP does not physically\n
                            locate within the region\n
\nCovariate:\n
    --cov-col       | -c    Header of covariates. If not provided, will use\n
                            all variables in the covariate file. By adding\n
                            @ in front of the string, any numbers within [\n
                            and ] will be parsed. E.g. @PC[1-3] will be\n
                            read as PC1,PC2,PC3. Discontinuous input are also\n
                            supported: @cov[1.3-5] will be parsed as \n
                            cov1,cov3,cov4,cov5\n
    --cov-factor            Header of categorical covariate(s). Dummy variable\n
                            will be automatically generated. Any items in\n
                            --cov-factor must also be found in --cov-col\n
                            Also accept continuous input (start with @).\n
    --cov-file      | -C    Covariate file. First column should be FID and \n
                            the second column should be IID. If --ignore-fid\n
                            is set, first column should be IID\n
\nP-value Thresholding:\n
    --bar-levels            Level of barchart to be plotted. When --fastscore\n
                            is set, PRSice will only calculate the PRS for \n
                            threshold within the bar level. Levels should be\n
                            comma separated without space\n
    --fastscore             Only calculate threshold stated in --bar-levels\n
    --no-full               By default, PRSice will include the full model, \n
                            i.e. p-value threshold = 1. Setting this flag will\n
                            disable that behaviour\n
    --interval      | -i    The step size of the threshold. Default: 0.00005 \n
    --lower         | -l    The starting p-value threshold. Default: 0.0001 \n
    --model                 Genetic model use for regression. The genetic\n
                            encoding is based on the base data where the\n
                            encoding represent number of the coding allele\n
                            Available models include:\n
                            add - Additive model, code as 0/1/2 (default)\n
                            dom - Dominant model, code as 0/1/1\n
                            rec - Recessive model, code as 0/0/1\n
                            het - Heterozygous only model, code as 0/1/0\n
    --no-regress            Do not perform the regression analysis and simply\n
                            output all PRS.\n
    --missing               Method to handle missing genotypes. By default, \n
                            final scores are averages of valid per-allele \n
                            scores with missing genotypes contribute an amount\n
                            proportional to imputed allele frequency. To throw\n
                            out missing observations instead (decreasing the\n
                            denominator in the final average when this happens),\n
                            use the 'no_mean_imputation' modifier. Alternatively,\n
                            you can use the 'center' modifier to shift all scores\n
                            to mean zero. \n
    --score                 Method to calculate the polygenic score.\n
                            Available methods include:\n
                            avg - Take the average effect size (default)\n
                            std - Standardize the effect size \n
                            sum - Direct summation of the effect size \n
    --upper         | -u    The final p-value threshold. Default: 0.5 \n
\nPRSet:\n
    --bed           | -B    Bed file containing the selected regions.\n
                            Name of bed file will be used as the region\n
                            identifier. WARNING: Bed file is 0-based\n
    --feature               Feature(s) to be included from the gtf file.\n
                            Default: exon,CDS,gene,protein_coding.\n
    --gtf           | -g    GTF file containing gene boundaries. Required\n
                            when --msigdb is used\n
    --msigdb        | -m    MSIGDB file containing the pathway information.\n
                            Require the gtf file\n
    --snp-set               Provide a SNP set file containing a single snp set.\n
                            Name of SNP set file will be used as the region\n
                            identifier. This file should contain only one column.\n
    --snp-sets              Provide a SNP set file containing multiple snp sets.\n
                            Each row represent a single SNP set with the first\n
                            column containing name of the SNP set.\n    
\nPRSlice:\n
    --prslice               Perform PRSlice where the whole genome is first cut\n
                            into bin size specified by this option. PRSice will\n
                            then be performed on each bin. Bins are then sorted\n
                            according to the their R2. PRSice is then performed\n
                            again to find the best bin combination.\n
                            This cannot be performed together with PRSet\n
                            (Currently not implemented)\n
\nPlotting:\n
    --bar-col-high          Colour of the most predicting threshold\n
                            Default: firebrick\n
    --bar-col-lower         Colour of the poorest predicting threshold\n
                            Default: dodgerblue\n
    --bar-col-p             Change the colour of bar to p-value threshold\n
                            instead of the association with phenotype\n
    --bar-palatte           Colour palatte to be used for bar plotting when\n
                            --bar_col_p is set. Default: YlOrRd\n
    --multi-plot            Plot the top N phenotype / gene set in a\n
                            summary plot\n 
    --plot                  When set, will only perform plotting.\n
    --plot-set              Define the gene set to be plot. Default: Base\n
    --quantile      | -q    Number of quantiles to plot. No quantile plot\n
                            will be generated when this is not provided.\n
    --quant-break           Quantile groupings for plotting the strata plot\n
    --quant-extract | -e    File containing sample ID to be plot on a separated\n
                            quantile e.g. extra quantile containing only \n
                            schizophrenia samples. Must contain IID. Should\n
                            contain FID if --ignore-fid isn't set.\n
    --quant-ref             Reference quantile for quantile plot\n
    --scatter-r2            y-axis of the high resolution scatter plot should be R2\n

\nMisc:\n
    --all-score             Output PRS for ALL threshold. WARNING: This\n
                            will generate a huge file\n
    --non-cumulate          Calculate non-cumulative PRS. PRS will be reset\n
                            to 0 for each new P-value threshold instead of\n
                            adding up\n
    --exclude               File contains SNPs to be excluded from the\n
                            analysis\n
    --extract               File contains SNPs to be included in the \n
                            analysis\n
    --ignore-fid            Ignore FID for all input. When this is set,\n
                            first column of all file will be assume to\n
                            be IID instead of FID\n
    --keep-ambig            Keep ambiguous SNPs. Only use this option\n
                            if you are certain that the base and target\n
                            has the same A1 and A2 alleles\n
    --logit-perm            When performing permutation, still use logistic\n
                            regression instead of linear regression. This\n
                            will substantially slow down PRSice\n
    --no-install            Forbid PRSice from automatically installing\n
                            the required packages (e.g. ggplot2)\n
    --out           | -o    Prefix for all file output\n
    --perm                  Number of permutation to perform. This swill\n
                            generate the empirical p-value. Recommend to\n
                            use value larger than 10,000\n
    --print-snp             Print all SNPs used to construct the best PRS\n
    --seed          | -s    Seed used for permutation. If not provided,\n
                            system time will be used as seed. When same\n
                            seed and same input is provided, same result\n
                            can be generated\n
    --thread        | -n    Number of thread use\n
    --x-range               Range of SNPs to be excluded from the whole\n
                            analysis. It can either be a single bed file\n
                            or a comma seperated list of range. Range must\n
                            be in the format of chr:start-end or chr:coordinate\n
    --help          | -h    Display this help message\n"




# Library handling --------------------------------------------------------

if (!exists('startsWith', mode = 'function')) {
    startsWith <- function(x, prefix) {
        return(substring(x, 1, nchar(prefix)) == prefix)
    }
}


libraries <-
    c("ggplot2",
      "data.table",
      "optparse",
      "methods",
      "tools",
      "grDevices",
      "RColorBrewer")
found.library.dir <- FALSE
argv <- commandArgs(trailingOnly = TRUE)
dir.arg.idx <- grep("--dir",argv)
no.install <- length(grep("--no-install", argv))>0
if (length(dir.arg.idx) != 0) {
    dir.arg.idx <- dir.arg.idx + 1
    found.library.dir <- TRUE
}

# INSTALL_PACKAGE: Functions for automatically install all required packages
InstalledPackage <- function(package) {
    available <- suppressMessages(suppressWarnings(
        sapply(
            package,
            require,
            quietly = TRUE,
            character.only = TRUE,
            warn.conflicts = FALSE
        )
    ))
    missing <- package[!available]
    if (length(missing) > 0)
        return(FALSE)
    return(TRUE)
}

CRANChoosen <- function()
{
    return(getOption("repos")["CRAN"] != "@CRAN@")
}

UsePackage <- function(package, dir, no.install)
{
    if (!InstalledPackage(package))
    {
        dir.create(file.path(dir, "lib"), showWarnings = FALSE)
        .libPaths(c(.libPaths(), paste(dir, "/lib", sep = "")))
        if (!InstalledPackage(package) & !no.install) {
            if (is.na(dir)) {
                writeLines("WARNING: dir not provided, cannot install the required packages")
                return(FALSE)
                
            } else{
                writeLines(paste(
                    "Trying to install ",
                    package,
                    " in ",
                    dir,
                    "/lib",
                    sep = ""
                ))
            }
            suppressMessages(suppressWarnings(
                install.packages(
                    package,
                    lib = paste(dir, "/lib", sep = ""),
                    repos = "http://cran.rstudio.com/",
                    quiet = T
                )
            ))
        }
        if (!InstalledPackage(package))
            return(FALSE)
    }
    return(TRUE)
}

use.data.table <- T
use.ggplot <- T #cerr
for (library in libraries)
{
    package.directory <- "."
    if (found.library.dir) {
        package.directory <- argv[dir.arg.idx]
    }
    if (!UsePackage(library, package.directory, no.install))
    {
        if (library == "data.table") {
            use.data.table <- F
            writeLines("Cannot install data.table, will fall back and use read.table instead")
            writeLines("Note: It will be slower when reading large files")
        } else if (library == "ggplot2") {
            use.ggplot <- F
            writeLines("Cannot install ggplot2, will fall back and native plotting devices")
            writeLines("Note: The legends will be uglier")
        } else{
            stop("Error: ", library, " cannot be load nor install!")
        }
    }
    
}

# Command line arguments --------------------------------------------------

# We don't type the help message here. Will just directly use the usage information from c++
# See the Help Messages section for more information
option_list <- list(
  # Base file
  make_option(c("--A1"), type = "character"),
  make_option(c("--A2"), type = "character"),
  make_option(c("-b", "--base"), type = "character"),
  make_option(c("--beta"), action = "store_true"),
  make_option(c("--bp"), type = "character"),
  make_option(c("--chr"), type = "character"),
  make_option(c("--index"), action = "store_true"),
  make_option(c("--info-base"), type = "character", dest = "info_base"),
  make_option(c("--maf-base"), type = "character", dest="maf_base"),
  make_option(c("--no-default"), action = "store_true", dest="no_default"),
  make_option(c("-p", "--pvalue"), type = "character"),
  make_option(c("--se"), type = "character"),
  make_option(c("--snp"), type = "character"),
  make_option(c("--stat"), type = "character"),
  # Target file
  make_option(c("--binary-target"), type = "character", dest = "binary_target"),
  make_option(c("--geno"), type = "numeric"),
  make_option(c("--info"), type = "numeric"),
  make_option(c("--keep"), type = "character"),
  make_option(c("--maf"), type = "numeric"),
  make_option(c("--nonfounders"), action = "store_true", dest = "nonfounders"),
  make_option(c("--pheno-col"), type = "character", dest = "pheno_col"),
  make_option(c("-f", "--pheno-file"), type = "character", dest = "pheno_file"),
  make_option(c("-k", "--prevalence"), type = "numeric"),
  make_option(c("--remove"), type = "character"),
  make_option(c("-t", "--target"), type = "character"),
  make_option(c("--target-list"), type = "character", dest="target_list"),
  make_option(c("--type"), type = "character"),
  # Dosage
  make_option(c("--allow-inter"), action = "store_true", dest="allow_inter"),
  make_option(c("--hard-thres"), type = "numeric", dest="hard_thres"),
  make_option(c("--hard"), action = "store_true"),
  # Clumping
  make_option(c("--clump-kb"), type = "character", dest = "clump_kb"),
  make_option(c("--clump-r2"), type = "numeric", dest = "clump_r2"),
  make_option(c("--clump-p"), type = "numeric", dest = "clump_p"),
  make_option(c("-L", "--ld"), type = "character"),
  make_option(c("--ld-list"), type = "character", dest="ld_list"),
  make_option(c("--ld-geno"), type = "numeric", dest="ld_geno"),
  make_option(c("--ld-info"), type = "numeric", dest="ld_info"),
  make_option(c("--ld-hard-thres"), type = "numeric", dest="ld_hard_thres"),
  make_option(c("--ld-keep"), type = "character", dest="ld_keep"),
  make_option(c("--ld-maf"), type = "numeric", dest="ld_maf"),
  make_option(c("--ld-remove"), type = "character", dest="ld_remove"),
  make_option(c("--ld-type"), type = "character", dest="ld_type"),
  make_option(c("--no-clump"), action = "store_true", dest = "no_clump"),
  make_option(c("--proxy"), type = "numeric"),
  # Covariates
  make_option(c("-c", "--cov-col"), type = "character", dest = "cov_col"),
  make_option(c("-C", "--cov-file"), type = "character", dest = "cov_file"),
  make_option(c("--cov-factor"), type = "character", dest = "cov_factor"),
  # P-thresholding
  make_option(
    c("--bar-levels"),
    type = "character",
    dest = "bar_levels"
  ),
  make_option(c("--fastscore"), action = "store_true"),
  make_option(c("--no-full"), action = "store_true"),
  make_option(c("-i", "--interval"), type = "numeric"),
  make_option(c("-l", "--lower"), type = "numeric"),
  make_option(c("--model"), type = "character"),
  make_option(c("--missing"), type = "character"),
  make_option(c("--no-regress"), action = "store_true", dest = "no_regress"),
  make_option(c("--score"), type = "character"),
  make_option(c("-u", "--upper"), type = "numeric"),
  # PRSet
  make_option(c("-B", "--bed"), type = "character"),
  make_option(c("--feature"), type = "character"),
  make_option(c("-g", "--gtf"), type = "character"),
  make_option(c("-m", "--msigdb"), type = "character"),
  make_option(c("--set-perm"), type = "numeric",dest="set_perm"),
  make_option(c("--wind-5"), type = "character", dest="wind_5"),
  make_option(c("--wind-3"), type = "character", dest="wind_3"),
  make_option(c("--snp-set"), type = "character", dest="snp_set"),
  make_option(c("--snp-sets"), type = "character", dest="snp_sets"),
  # PRSlice 
  make_option(c("--prslice"), type = "numeric"),
  # Misc
  make_option(c("--all-score"), action = "store_true", dest="all_score"),
  make_option(c("--exclude"), type = "character"),
  make_option(c("--extract"), type = "character"),
  make_option(c("--ignore-fid"), action = "store_true", dest = "ignore_fid"),
  make_option(c("--logit-perm"), action = "store_true", dest = "logit_perm"),
  make_option(c("--keep-ambig"), action = "store_true", dest = "keep_ambig"),
  make_option(c("-o", "--out"), type = "character", default = "PRSice"),
  make_option(c("--perm"), type = "numeric"),
  make_option(c("-s", "--seed"), type = "numeric"),
  make_option(c("--print-snp"), action = "store_true", dest = "print_snp"),
  make_option(c("--non-cumulate"), action = "store_true", dest = "non_cumulate"),
  make_option(c("--pearson"), action = "store_true"),
  make_option(c("-n", "--thread"), type = "numeric"),
  make_option(c("--x-range"), type = "character", dest="x_range"),
  #R Specified options
  make_option(c("--plot"), action = "store_true"),
  make_option(c("--quantile", "-q"), type = "numeric"),
  make_option(c("--quant-break"), type="character", dest="quant_break"),
  make_option(c("--multi-plot"), type = "numeric", dest = "multi_plot"),
  make_option(c("--plot-set"),
              type = "character",
              dest = "plot_set",
              default = "Base"),
  make_option(c("--quant-pheno"), action = "store_true", dest = "quant_pheno"),
  make_option(c("--quant-extract", "-e"), type = "character", dest = "quant_extract"),
  make_option("--quant-ref", type = "numeric", dest = "quant_ref"),
  make_option("--scatter-r2",
              action = "store_true",
              default = F,
              dest = "scatter_r2"),
  make_option("--bar-col-p",
              action = "store_true",
              default = F,
              dest = "bar_col_p"),
  make_option("--bar-col-low",
              type = "character",
              default = "dodgerblue",
              dest = "bar_col_low"),
  make_option("--bar-col-high",
              type = "character",
              default = "firebrick",
              dest = "bar_col_high"),
  make_option("--bar-palatte",
              type = "character",
              default = "YlOrRd",
              dest = "bar_palatte"),
  make_option("--prsice", type = "character"),
  make_option("--dir", type = "character")
)


capture <- commandArgs(trailingOnly = TRUE)
help <- (sum(c("--help", "-h") %in% capture) >= 1)
has_c <- (sum(c("--prsice") %in% capture) >= 1)
if (help) {
    cat(help_message)
    quit()
}
argv <- parse_args(OptionParser(option_list = option_list))

not_cpp <- c(
    "help",
    "plot",
    "quantile",
    "quant-extract",
    "intermediate",
    "quant-ref",
    "quant-pheno",
    "quant-break",
    "scatter-r2",
    "bar-col-p",
    "bar-col-low",
    "bar-col-high",
    "bar-palatte",
    "prsice",
    "multi-plot",
    "plot-set",
    "dir",
    "no-install"
)

if (is.null(argv$cov_col) && !is.null(argv$cov_header))
{
    argv$cov_col = argv$cov_header
}

# Check help messages --------------------------------------------------

provided <- function(name, argv) {
    return(name %in% names(argv))
}

get_os <- function(){
    sysinf <- Sys.info()
    if (!is.null(sysinf)){
        os <- sysinf['sysname']
        if (os == 'Darwin')
            os <- "osx"
    } else { ## mystery machine
        os <- .Platform$OS.type
        if (grepl("^darwin", R.version$os))
            os <- "osx"
        if (grepl("linux-gnu", R.version$os))
            os <- "linux"
    }
    tolower(os)
}

# CALL_PRSICE: Call the cpp PRSice if required
# To ensure the excutable is set correctly
# For window, we might not be able to start an executable by simply adding ./
# therefore Window people will need to be careful with their parameter input
os <- get_os()
if (provided("prsice", argv)) {
    if (!startsWith(argv$prsice, "/") &&
        !startsWith(argv$prsice, ".") && os!="windows") {
        argv$prsice = paste("./", argv$prsice, sep = "")
    }
}

# Running PRSice ----------------------------------------------------------

# We don't bother to check if the input is correct, the parameter should be checked by the c++ program
add_command <- function(input) {
    if (length(input) == 1) {
        if (is.na(input)) {
            return(NA)
        } else{
            return(input)
        }
    } else{
        return(paste(input, collapse = ","))
    }
}
command <- ""
argv_c <- argv

names(argv_c) <- gsub("_", "-", names(argv))
flags <-
    c(
        "all-score",
        "allow-inter",
        "beta",
        "fastscore",
        "ignore-fid",
        "index",
        "keep-ambig",
        "logit-perm",
        "no-clump",
        "no-default",
        "no-full",
        "no-mt",
        "no-regress",
        "no-x",
        "no-xy",
        "no-y",
        "non-cumulate",
        "print-snp"
    )

if (!provided("plot", argv)) {
    for (i in names(argv_c)) {
        # only need special processing for flags and specific inputs
        if (i %in% flags) {
            if (argv_c[[i]])
                command = paste(command, " --", i, sep = "")
        } else if (i %in% not_cpp) {
            # ignore
        } else{
            temp = add_command(argv_c[[i]])
            if (!is.na(temp)) {
                command = paste(command, " --", i, " ", temp, sep = "")
            }
        }
    }
    if (nchar(command) == 0) {
        cat(help_message)
        quit()
    }
    if (provided("prsice", argv_c)) {
        ret <- system2(argv_c$prsice, command)
        if (ret != 0 || provided(help, argv)) {
            stop()
            
        }
    } else{
        stop("Cannot run PRSice without the PRSice binary file")
    }
}

# Helper functions --------------------------------------------------------
max_length <- function(x) {
    info <- strsplit(as.character(x), split = "\n")[[1]]
    max(sapply(info, nchar))
}
str_wrap <- function(x) {
    lapply(strwrap(x, width = 25, simplify = FALSE), paste, collapse = "\n")
}
shorten_label <- function(x) {
    lab <-
        paste(strsplit(paste(
            strsplit(as.character(x), split = "\\.")[[1]], collapse = " "
        ), split = "_")[[1]], collapse = " ")
    return(str_wrap(lab)[[1]])
}


get_quantile <- function(x, num.quant, quant.ref){
    quant <- as.numeric(cut(x,
                            breaks = unique(quantile(
                                x, probs = seq(0, 1, 1 / num.quant)
                            )),
                            include.lowest = T))
    if(is.na(quant.ref) | is.null(quant.ref)){
        quant.ref <- ceiling(num.quant / 2)
    }
    quant <- factor(quant, levels = c(quant.ref, seq(min(quant), max(quant), 1)[-quant.ref]))
    return(quant)
}

set_uneven_quant <- function(quant.cutoff, ref.cutoff, num.quant, prs, quant.index){

    quant <- get_quantile(prs, num.quant, 1)
    quant <- factor(quant,1:num.quant)
    quant.cut <- sort(as.numeric(strsplit(quant.cutoff, split=",")[[1]]))
    if((!is.null(ref.cutoff) & sum(ref.cutoff == quant.cut)==0) | num.quant < max(quant.cut)){
        stop(
            "Invalid quant-break. quant-break must be smaller than total number of quantiles and quant-ref must be one of the quant-break"
        )
    }
    prev.name <- 0
    ref.level <- NULL
    for(i in 1:length(quant.cut)){
        up.bound <- max(which(suppressWarnings(as.numeric(levels(quant))) <= quant.cut[i]))
        cur.name <- levels(quant)[up.bound]
        name.level <- paste0("(",prev.name, ",", cur.name, "]")
        
        if(prev.name==0){
            name.level <- paste0("[",prev.name, ",", cur.name, "]")
            
        }
        if(!is.null(ref.cutoff)) {
            if (quant.cut[i] == ref.cutoff) {
                ref.level <- name.level
                
            }
        }
        range <- i:up.bound
        levels(quant)[range] <- rep(name.level, length(range))
        prev.name <- cur.name
    }
    if(is.null(ref.level)){
        writeLines("=======================================")
        writeLines("Warning: Cannot find required reference level, will use the middle level as reference")
        writeLines("=======================================")
        ref.level <- levels(quant)[ceiling(length(quant.cut)/2)]
    }
    ref.index <- which(levels(quant)==ref.level)
    quant.index <- c(ref.index,c(1:length(levels(quant)))[-ref.index])
    quant<- relevel(quant, ref=ref.level)
    return(list(quant, quant.index))
}
# Determine Default -------------------------------------------------------
# First, determine the bar levels
if(!provided("bar_levels", argv)){
    if(!provided("msigdb", argv) & !provided("gtf", argv) & !provided("bed", argv)) {
        argv$bar_levels <- paste(0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, sep=",")
        if (!provided("no_full", argv)) {
            argv$bar_levels <- paste(argv$bar_levels, 1, sep=",")
        }
    } else if (!provided("fastscore", argv) &
               !provided("lower", argv) &
               !provided("upper", argv) & !provided("interval", argv)) {
        # This is prset, so by default, we don't do all threshold
        # unless user use some of the parameter related to the thresholding
        argv$bar_levels <- "1"
    }
}
# Next, we need to determine if we are doing binary target
# This is only required when user does not provide --binary-target
 
if(!provided("binary_target", argv)){
    # Now we want to check if base is beta
    base_beta <- F
    if(provided("beta", argv)){
        # Base is beta
        base_beta <- T
    }else{
        if(!provided("base", argv)){
            writeLines("Warning: Without base file, we cannot determine if the summary statistic is beta or OR, which is used to determine if the target phenotype is binary or not. We will now proceed assuming the target phenotype is binary. If that is incorrect, please use --binary-target to specify the correct target phenotype type, or you can provide the base file")
        }
        else{
            zz <- gzfile(argv$base)
            base_header <- readLines(zz,n=1)
            close(zz)
            or <- length(grep("or",base_header,ignore.case=TRUE))==1
            beta <- length(grep("beta",base_header,ignore.case=TRUE))==1
            if(or & beta){
                stop("Both OR and BETA detected. Cannot determine which one should be used. Please use --beta or --binary-target")
            }
            if(beta){
                base_beta <- T    
            }
            if(!or & !beta){
                stop("Do not detect either BETA or OR. Please ensure your base input is correct")
            }
        }
    }
    # Now we know if the base is beta or not, we can determine the target binary status
    if(!provided("pheno_col", argv)){
        if(base_beta){
            argv$binary_target <- "F"
        }else{
            argv$binary_target <- "T"
        }
    }else if(length(strsplit(argv$pheno_col, split=",")[[1]])==1){
        if(base_beta){
            argv$binary_target <- "F"
        }else{
            argv$binary_target <- "T"
        }    
    }
}

# Sanity check for binary-target
if(provided("pheno_col", argv) & provided("binary_target", argv)){
    pheno_length <- length(strsplit(argv$pheno_col, split=",")[[1]])
    binary_length <- length(strsplit(argv$binary_target, split=",")[[1]])
    if(pheno_length==0 & binary_length==1){
        # This is ok
    }else if(pheno_length!=binary_length){
        stop("Error: Number of target phenotypes doesn't match information of binary target! You must indicate whether the phenotype is binary using --binary-target\n")        
    }
}
    

# Plottings ---------------------------------------------------------------

# Standard Theme for all plots
theme_sam <- NULL
if(use.ggplot){
  theme_sam <- theme_bw()+theme(axis.title=element_text(face="bold", size=18),
                              axis.text=element_text(size=14),
                              legend.title=element_text(face="bold", size=18),
                              legend.text=element_text(size=14),
                              axis.text.x=element_text(angle=45, hjust=1),
                              panel.grid = element_blank(),
                              panel.border = element_blank(), 
                              axis.line = element_line()
                              )
}

# Quantile Plots----------------------------------------------------------

call_quantile <-
    function(pheno.merge,
             prefix,
             num_quant,
             quant.index,
             pheno.as.quant,
             use.residual,
             use.ggplot,
             binary,
             extract,
             uneven) {
        
    if (!pheno.as.quant) {
        family <- gaussian
        if (binary) {
            if (!use.residual) {
                family <- binomial
            }
        }
        reg <-
            summary(glm(Pheno ~ quantile, family, data = pheno.merge))
        coef.quantiles <- (reg$coefficients[1:num_quant, 1])
        ci <- (1.96 * reg$coefficients[1:num_quant, 2])
        
        ci.quantiles.u <-
            coef.quantiles + ci
        ci.quantiles.l <-
            coef.quantiles - ci
        if (binary & !use.residual) {
            ci.quantiles.u <- exp(ci.quantiles.u)
            ci.quantiles.l <- exp(ci.quantiles.l)
            coef.quantiles <- exp(coef.quantiles)
        }
        
        coef.quantiles[1] <-
            ifelse(binary & !use.residual, 1, 0)
        ci.quantiles.u[1] <-
            ifelse(binary & !use.residual, 1, 0)
        ci.quantiles.l[1] <-
            ifelse(binary & !use.residual, 1, 0)
        quantiles.for.table <- factor(levels(pheno.merge$quantile), levels(pheno.merge$quantile))
        quantiles.df <-
            data.frame(
                Coef = coef.quantiles,
                CI.U = ci.quantiles.u,
                CI.L = ci.quantiles.l,
                DEC = quantiles.for.table
            )
        quantiles.df$Group = 0
        if (!is.null(extract)) {
            # Last element should be the cases
            quantiles.df$Group[nrow(quantiles.df)] <- 1
        }
        quantiles.df$Group <-
            factor(quantiles.df$Group, levels = c(0, 1))
        quantiles.df <- quantiles.df[order(quant.index),]
        quantiles.df$DEC <- factor(quantiles.df$DEC, levels=levels(quantiles.df$DEC)[order(quant.index)])
        row.names(quantiles.df) <- quantiles.df$DEC
        quantiles.df <- cbind(Quantile = rownames(quantiles.df), quantiles.df) 
        quant.out <- quantiles.df[,-5]
        sample.size <- as.data.frame(table(pheno.merge$quantile))
        colnames(sample.size) <- c("Quantile", "N")
        quant.out$Order <- 1:nrow(quant.out)
        quant.out <- merge(quant.out,sample.size )
        if(binary & !use.residual){
            colnames(quant.out)[2] <- "OR"
        }
        
        quant.out <-
            quant.out[order(quant.out$Order), !colnames(quant.out) %in% "Order"]
        write.table(
            quant.out,
            paste(prefix, "_QUANTILES_", Sys.Date(), ".txt", sep = ""),
            sep = "\t",
            quote = F,
            row.names = F
        )
        if (use.ggplot) {
            plot.quant(quantiles.df,
                       num_quant,
                       binary,
                       extract,
                       prefix,
                       use.residual,
                       uneven)
        } else{
            plot.quant.no.g(quantiles.df,
                            num_quant,
                            binary,
                            extract,
                            prefix,
                            use.residual,
                            uneven)
        }
    } else{
        # TODO: Maybe also change this to regression? Though might be problematic if we have binary pheno without cov
        pheno.sum <-
            data.frame(
                mean = numeric(num_quant),
                quantile = factor(levels(pheno.merge$quantile)[order(quant.index)],levels=levels(pheno.merge$quantile)[order(quant.index)]),
                UCI = numeric(num_quant),
                LCI = numeric(num_quant)
            )
        for (i in 1:length(levels(pheno.sum$quantile))) {
            
            cur.prs <-
                pheno.merge$PRS[as.character(pheno.merge$quantile) %in% as.character(levels(pheno.sum$quantile))[i]]
            pheno.sum$mean[i] <- mean(cur.prs, na.rm = T)
            pheno.sum$UCI[i] <-
                pheno.sum$mean[i] + sd(cur.prs, na.rm = T)
            pheno.sum$LCI[i] <-
                pheno.sum$mean[i] - sd(cur.prs, na.rm = T)
        }
        pheno.sum$Group = 0
        if (!is.null(extract)) {
            pheno.sum$Group[num_quant] = 1
        }
        pheno.sum$Group <-
            factor(pheno.sum$Group, levels = c(0, 1))
        write.table(
            pheno.sum,
            paste(prefix, "_PHENO_QUANTILES_", Sys.Date(), ".txt", sep = ""),
            sep = "\t",
            quote = F,
            row.names = F
        )
        if (use.ggplot) {
            plot.pheno.quant(pheno.sum,
                             use.residual,
                             num_quant,
                             extract,
                             prefix,
                             uneven)
        } else{
            plot.pheno.quant.no.g(pheno.sum,
                                  use.residual,
                                  num_quant,
                                  extract,
                                  prefix,
                                  uneven)
        }
    }
    
    }


uneven_quantile_plot <- function(base.prs, pheno, prefix, argv, binary, use.ggplot, use.residual){
    binary <- as.logical(binary)
    writeLines("Plotting the quantile plot")
    extract <- NULL
    if (provided("quant_extract", argv)) {
        if (use.data.table) {
            extract <- fread(argv$quant_extract,
                             header = F,
                             data.table = F)
        } else{
            extract <- read.table(argv$quant_extract,
                                  header = F)
        }
        
    }
    num_quant <- argv$quantile
    pheno.merge <- merge(base.prs, pheno)
    pheno.as.quant <- provided("quant_pheno", argv)
    if (pheno.as.quant &&
        length(unique(pheno.merge$Pheno)) < num_quant) {
        writeLines(
            paste(
                "WARNING: There are only ",
                length(unique(pheno.merge$Pheno)),
                " unique Phenotype but asked for ",
                num_quant,
                " quantiles",
                sep = ""
            )
        )
        writeLines(paste("Will not generate the quantile plot for ", prefix))
        return()
    } else if (length(unique(pheno.merge$PRS)) < num_quant) {
        writeLines(
            paste(
                "WARNING: There are only ",
                length(unique(pheno.merge$PRS)),
                " unique PRS but asked for ",
                num_quant,
                " quantiles",
                sep = ""
            )
        )
        writeLines(paste("Will not generate the quantile plot for ", prefix))
        return()
    }
    quants <- NULL
    quant.index <- NULL
    if (!pheno.as.quant) {
        quant.info <- set_uneven_quant(argv$quant_break, argv$quant_ref, num_quant, pheno.merge$PRS, quant.index)
        quants <- quant.info[[1]]
        quant.index <- quant.info[[2]]
    } else{
        quant.info <- set_uneven_quant(argv$quant_break, argv$quant_ref, num_quant, pheno.merge$Pheno, quant.index)
        quants <- quant.info[[1]]
        quant.index <- quant.info[[2]]
    }
    
    num_quant <- length(levels(quants))
    if (!is.null(extract)) {
        extract_ID <- NULL
        best_ID <- NULL
        if (provided("ignore_fid", argv)) {
            extract_ID <- extract$V1
            best_ID <- pheno.merge$IID
        } else{
            extract_ID <- paste(extract$V1, extract$V2, sep = "_")
            best_ID <-
                paste(pheno.merge$FID, pheno.merge$IID, sep = "_")
        }
        levels(quants) <- c(levels(quants),"Extracted")
        quants[best_ID %in% extract_ID] <- "Extracted"
        num_quant <- num_quant + 1
    }
    pheno.merge$quantile <- quants
    call_quantile(pheno.merge,
                  prefix,
                  num_quant,
                  quant.index,
                  pheno.as.quant,
                  use.residual,
                  use.ggplot,
                  binary,
                  extract, TRUE)
}


quantile_plot <- function(base.prs, pheno, prefix, argv, binary, use.ggplot, use.residual){
    binary <- as.logical(binary)
    writeLines("Plotting the quantile plot")
    extract <- NULL
    if (provided("quant_extract", argv)) {
        if (use.data.table) {
            extract <- fread(argv$quant_extract,
                             header = F,
                             data.table = F)
        } else{
            extract <- read.table(argv$quant_extract,
                                  header = F)
        }
        
    }
    num_quant <- argv$quantile
    
    pheno.merge <- merge(base.prs, pheno)
    pheno.as.quant <- provided("quant_pheno", argv)
    if (pheno.as.quant &&
        length(unique(pheno.merge$Pheno)) < num_quant) {
        writeLines(
            paste(
                "WARNING: There are only ",
                length(unique(pheno.merge$Pheno)),
                " unique Phenotype but asked for ",
                num_quant,
                " quantiles",
                sep = ""
            )
        )
        writeLines(paste("Will not generate the quantile plot for ", prefix))
        return()
    } else if (length(unique(pheno.merge$PRS)) < num_quant) {
        writeLines(
            paste(
                "WARNING: There are only ",
                length(unique(pheno.merge$PRS)),
                " unique PRS but asked for ",
                num_quant,
                " quantiles",
                sep = ""
            )
        )
        writeLines(paste("Will not generate the quantile plot for ", prefix))
        return()
    }
    
    quant.ref <- ceiling(argv$quantile / 2)
    if (provided("quant_ref", argv)) {
        quant.ref <- argv$quant_ref
        if (quant.ref > argv$quantile) {
            quant.ref <- ceiling(argv$quantile / 2)
            writeLines("=======================================")
            writeLines(
                paste(
                    "WARNING: reference quantile",
                    quant.ref,
                    "is greater than number of quantiles",
                    argv$quantile,
                    "\n Using middle quantile by default"
                )
            )
            writeLines("=======================================")
        }
    }
    quants <- NULL
    if (!pheno.as.quant) {
        quants <- get_quantile(pheno.merge$PRS, num_quant, quant.ref)
    } else{
        quants <- get_quantile(pheno.merge$Pheno, num_quant, quant.ref)
    }
    num_quant <- length(levels(quants))
    if (!is.null(extract)) {
        extract_ID <- NULL
        best_ID <- NULL
        if (provided("ignore_fid", argv)) {
            extract_ID <- extract$V1
            best_ID <- pheno.merge$IID
        } else{
            extract_ID <- paste(extract$V1, extract$V2, sep = "_")
            best_ID <-
                paste(pheno.merge$FID, pheno.merge$IID, sep = "_")
        }
        levels(quants) <- c(levels(quants),"Extracted")
        quants[best_ID %in% extract_ID] <- "Extracted"
        num_quant <- num_quant + 1
    }
    quant.index <- c(quant.ref, c(1:num_quant)[-quant.ref])
    pheno.merge$quantile <- quants
    call_quantile(pheno.merge,
                 prefix,
                 num_quant,
                 quant.index,
                 pheno.as.quant,
                 use.residual,
                 use.ggplot,
                 binary,
                 extract, FALSE)
}

plot.pheno.quant.no.g <- function(pheno.sum, use_residual, num_quant, extract, prefix, uneven){
  png(paste(prefix, "_QUANTILES_PHENO_PLOT_", Sys.Date(), ".png", sep = ""),
      height=10, width=10, res=300, unit="in")
  par(pty="s", cex.lab=1.5, cex.axis=1.25, font.lab=2, mai=c(0.5,1.25,0.1,0.1))
  pheno.sum$color <- "#D55E00"
  xlab <- NULL
  name <- "Quantiles"
  if(uneven){
      name <- "Strata"
  }
  if(use_residual){
    xlab <-paste0(name, " for Residualized Phenotype")
  }else{
    xlab <-paste0(name, " for Phenotype")
  }
  
  
  if(!is.null(extract)){
    pheno.sum$color <- "#D55E00"
    pheno.sum$color[quantiles.df$Group==1] <- "#0072B2"
  }
  ylab <- "Mean PRS given phenotype in quantiles"
  if(uneven){
      ylab <- "Mean PRS given phenotype in strata"
  }
  with(pheno.sum, 
       plot(x=quantile, y=mean, 
            col=color, pch=19, 
            axes=F, cex=1.5,
            ann=F,
            ylim=c(min(LCI),max(UCI))
       ))
  box(bty='L', lwd=2)
  axis(2,las=2, lwd=2)
  axis(1, label=seq(1,num_quant,2), at=seq(1,num_quant,2),lwd=2)
  axis(1, label=seq(2,num_quant,2), at=seq(2,num_quant,2),lwd=2)
  with(pheno.sum, arrows(quantile,mean, quantile,LCI,length=0, col=color, lwd=1.5))
  with(pheno.sum, arrows(quantile,mean, quantile,UCI,length=0, col=color, lwd=1.5))
  title(ylab=ylab, line=4, cex.lab=1.5, font=2 )
  title(xlab=xlab, line=2.5, cex.lab=1.5, font=2 )
  g<-dev.off()
}

plot.pheno.quant <- function(pheno.sum, use_residual, num_quant, extract, prefix, uneven){
  quantiles.plot <-
    ggplot(pheno.sum, aes(
      x = quantile,
      y = mean,
      ymin = LCI,
      ymax = UCI
    ))+ 
    theme_sam
  if(uneven) {
      quantiles.plot <-
          quantiles.plot + ylab("Mean PRS given phenotype in strata")
  } else{
      quantiles.plot <-
          quantiles.plot + ylab("Mean PRS given phenotype in quantiles")
  }
  if(use_residual){
      if (uneven) {
          quantiles.plot <- quantiles.plot +
              xlab("Strata for Residualized Phenotype")
      } else{
          quantiles.plot <- quantiles.plot +
              xlab("Quantiles for Residualized Phenotype")
      }
  }else{
      if (uneven) {
          quantiles.plot <- quantiles.plot +
              xlab("Quantiles for Phenotype")
      } else{
          quantiles.plot <- quantiles.plot +
              xlab("Strat for Phenotype")
      }
  }
  
  if (is.null(extract)) {
    quantiles.plot <-
      quantiles.plot + geom_point(colour = "#D55E00", size = 4) +
      geom_pointrange(colour = "#D55E00", size = 0.9)
  } else{
    quantiles.plot <-
      quantiles.plot + geom_point(aes(color = Group), size = 4) +
      geom_pointrange(aes(color = Group), size = 0.9) +
      scale_colour_manual(values = c("#D55E00","#0072B2"))
  }
  ggsave(
    paste(prefix, "QUANTILES_PHENO_PLOT_", Sys.Date(),".png", sep = "_"),
    quantiles.plot,
    height=10,width=10
  )
}

plot.quant <- function(quantiles.df, num_quant, binary, extract, prefix, use_residual, uneven){
    quantiles.plot <-
        ggplot(quantiles.df, aes(
            x = DEC,
            y = Coef,
            ymin = CI.L,
            ymax = CI.U
        )) + 
        theme_sam
    if(uneven){
        quantiles.plot <- quantiles.plot+xlab("Strata for Polygenic Score")
    } else{
        quantiles.plot <-
            quantiles.plot + xlab("Quantiles for Polygenic Score")
    }
    if (binary){
        if(!use_residual) {
            quantiles.plot <-
                quantiles.plot + ylab("Odds Ratio for Score on Phenotype")
        }else if(use_residual){
            if(uneven){
                quantiles.plot <- quantiles.plot + ylab("Change in residualized\nPhenotype given score in strata")
            }else{
                quantiles.plot <- quantiles.plot + ylab("Change in residualized\nPhenotype given score in quantiles")
            }
        }
    } else if(use_residual){
        if (uneven) {
            quantiles.plot <-
                quantiles.plot + ylab("Change in residualized\nPhenotype given score in strata")
        } else{
            quantiles.plot <-
                quantiles.plot + ylab("Change in residualized\nPhenotype given score in quantiles")
        }
    }else{
        if (uneven) {
            quantiles.plot <- quantiles.plot +
                ylab("Change in Phenotype \ngiven score in strata")
        } else{
            quantiles.plot <- quantiles.plot +
                ylab("Change in Phenotype \ngiven score in quantiles")
        }
    }
    if (is.null(extract)) {
        quantiles.plot <-
            quantiles.plot + geom_point(colour = "royalblue2", size = 4) +
            geom_pointrange(colour = "royalblue2", size = 0.9)
    } else{
        quantiles.plot <-
            quantiles.plot + geom_point(aes(color = Group), size = 4) +
            geom_pointrange(aes(color = Group), size = 0.9) +
            scale_colour_manual(values = c("#0072B2", "#D55E00"))
    }
    ggsave(
        paste(prefix, "_QUANTILES_PLOT_", Sys.Date(),".png", sep = ""),
        quantiles.plot,
        height=10, width=10
    )
}

plot.quant.no.g <- function(quantiles.df, num_quant, binary, extract, prefix, use_residual, uneven){
  png(paste(prefix, "_QUANTILES_PLOT_", Sys.Date(), ".png", sep = ""),
      height=10, width=10, res=300, unit="in")
  par(pty="s", cex.lab=1.5, cex.axis=1.25, font.lab=2, mai=c(0.5,1.25,0.1,0.1))
  quantiles.df$color <- "royalblue2"
  if(!is.null(extract)){
    quantiles.df$color <- "#0072B2"
    quantiles.df$color[quantiles.df$Group==1] <- "#D55E00"
  }
  ylab <- NULL
  if (binary){
      if(!use_residual) {
          quantiles.plot <-
              quantiles.plot + ylab("Odds Ratio for Score on Phenotype")
      }else if(use_residual){
          if(uneven) {
              quantiles.plot <-
                  quantiles.plot + ylab("Change in residualized\nPhenotype given score in strata")
          } else{
              quantiles.plot <-
                  quantiles.plot + ylab("Change in residualized\nPhenotype given score in quantiles")
          }
      }
  } else if(use_residual){
      if (uneven) {
          quantiles.plot <-
              quantiles.plot + ylab("Change in residualized\nPhenotype given score in strata")
      } else{
          quantiles.plot <-
              quantiles.plot + ylab("Change in residualized\nPhenotype given score in quantiles")
      }
  } else{
      if(uneven) {
          quantiles.plot <- quantiles.plot +
              ylab("Change in Phenotype \ngiven score in strata")
      } else{
          quantiles.plot <- quantiles.plot +
              ylab("Change in Phenotype \ngiven score in quantiles")
      }
  }
  xlab <- "Quantiles for Polygenic Score"
  if(uneven){
      xlab <- "Strata for Polygenic Score"
  }
  with(quantiles.df, 
       plot(x=DEC, y=Coef, 
            col=color, pch=19, 
            axes=F, cex=1.5, ann=F,
            ylim=c(min(CI.L),max(CI.U))
            ))

  axis(2,las=2,lwd=2)
  box(bty='L', lwd=2)
  axis(1, label=seq(1,num_quant,2), at=seq(1,num_quant,2),lwd=2)
  axis(1, label=seq(2,num_quant,2), at=seq(2,num_quant,2),lwd=2)
  with(quantiles.df, arrows(DEC,Coef, DEC,CI.L,length=0, col=color, lwd=1.5))
  with(quantiles.df, arrows(DEC,Coef, DEC,CI.U,length=0, col=color, lwd=1.5))
  title(ylab=ylab, line=4, cex.lab=1.5, font=2 )
  title(xlab=xlab, line=2.5, cex.lab=1.5, font=2 )
  g<-dev.off()
}


# High Resolution Plot ----------------------------------------------------


high_res_plot <- function(PRS, prefix, argv, use.ggplot) {
    # we will always include the best threshold
    writeLines("Plotting the high resolution plot")
    barchart.levels <-
        c(strsplit(argv$bar_levels, split = ",")[[1]], PRS$Threshold[which.max(PRS$R2)])
    barchart.levels <-
        as.numeric(as.character(sort(
            unique(barchart.levels), decreasing = F
        )))
    # As the C++ program will skip thresholds, we need to artificially add the correct threshold information
    PRS.ori <- PRS
    threshold <- as.numeric(as.character(PRS.ori$Threshold))
    for (i in 1:length(barchart.levels)) {
        if (sum(barchart.levels[i] - threshold > 0) > 0) {
            target <- max(threshold[barchart.levels[i] - threshold >= 0])
            temp <- PRS.ori[threshold == target,]
            temp$Threshold <- barchart.levels[i]
            PRS <- rbind(PRS, temp)
            
        } else{
            target <-
                (threshold[which(abs(threshold - barchart.levels[i]) == min(abs(threshold - barchart.levels[i])))])
            temp <- PRS.ori[threshold == target,]
            temp$Threshold <- barchart.levels[i]
            PRS <- rbind(PRS, temp)
            
        }
    }
    PRS = unique(PRS)
    # Need to also plot the barchart level stuff with green
    if(use.ggplot){
      plot.high.res(argv, PRS, prefix, barchart.levels)
    }else{
      plot.high.res.no.g(argv, PRS, prefix, barchart.levels)
    }
}

plot.high.res.no.g <- function(argv, PRS, prefix, barchart.levels){
  png(paste(prefix, "_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""),
      height=10, width=10, res=300, unit="in")
  par(pty="s", cex.lab=1.5, cex.axis=1.25, font.lab=2, mai=c(0.5,1.25,0.1,0.1))
  xlab <- expression(italic(P) - value ~ threshold ~ (italic(P)[T]))
  ylab <- expression(paste("PRS model fit:  ", R ^ 2, sep = " "))
  if (argv$scatter_r2) {
    with(PRS[order(PRS$Threshold),], 
         plot(x=Threshold, y=R2, 
              pch=20, 
              axes=FALSE, ann=F))
    
    with(subset(PRS[order(PRS$Threshold),], Threshold%in%barchart.levels ), 
         lines(x=Threshold, y=R2, col="green"))
  }else{
    ylab <- bquote(PRS ~ model ~ fit: ~ italic(P) - value ~ (-log[10]))
    with(PRS[order(PRS$Threshold),], 
         plot(x=Threshold, y=-log10(P),  
              pch=20, 
              axes=FALSE, ann=F))
    
    with(subset(PRS[order(PRS$Threshold),], Threshold%in%barchart.levels ), 
         lines(x=Threshold, y=-log10(P), col="green"))
  }
  box(bty='L', lwd=2)
  axis(2, las=2)
  axis(1)
  
  title(ylab=ylab, line=4, cex.lab=1.5, font=2 )
  title(xlab=xlab, line=2.5, cex.lab=1.5, font=2 )
  g<-dev.off()
}

plot.high.res <- function(argv, PRS, prefix, barchart.levels){
  ggfig.points <- ggplot(data = PRS, aes(x = Threshold)) +
    xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
    theme_sam
  if (argv$scatter_r2) {
    ggfig.points <-
      ggfig.points + geom_point(aes(y = R2)) + geom_line(aes(y = R2), colour = "green",
                                                         data = PRS[with(PRS, Threshold %in% barchart.levels) ,]) +
      ylab(expression(paste("PRS model fit:  ", R ^ 2, sep = " ")))
  } else{
    ggfig.points <-
      ggfig.points + geom_point(aes(y = -log10(P))) + geom_line(aes(y = -log10(P)), colour = "green",
                                                                data = PRS[with(PRS, Threshold %in% barchart.levels) ,]) +
      ylab(bquote(PRS ~ model ~ fit: ~ italic(P) - value ~ (-log[10])))
    
  }
  ggsave(
    paste(prefix, "_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""),
    ggfig.points,
    height=10, width=10
  )
}


# Plot bar plot -----------------------------------------------------------


bar_plot <- function(PRS, prefix, argv, use.ggplot) {
    writeLines("Plotting Bar Plot")
    barchart.levels <-
        c(strsplit(argv$bar_levels, split = ",")[[1]], PRS$Threshold[which.max(PRS$R2)])
    barchart.levels <-
        as.numeric(as.character(sort(
            unique(barchart.levels), decreasing = F
        )))
    threshold <- as.numeric(as.character(PRS$Threshold))
    PRS.ori = PRS
    threshold <- as.numeric(as.character(PRS.ori$Threshold))
    for (i in 1:length(barchart.levels)) {
        if (sum(barchart.levels[i] - threshold > 0) > 0) {
            target <- max(threshold[barchart.levels[i] - threshold >= 0])
            temp <- PRS.ori[threshold == target,]
            temp$Threshold <- barchart.levels[i]
            PRS = rbind(PRS, temp)
            
        } else{
            target <-
                (threshold[which(abs(threshold - barchart.levels[i]) == min(abs(threshold - barchart.levels[i])))])
            temp <- PRS.ori[threshold == target,]
            temp$Threshold <- barchart.levels[i]
            PRS = rbind(PRS, temp)
            
        }
    }
    PRS <- unique(PRS[order(PRS$Threshold),])
    # As the C++ program will skip thresholds, we need to artificially add the correct threshold information
    output <- PRS[PRS$Threshold %in% barchart.levels, ]
    output$print.p[round(output$P, digits = 3) != 0] <-
        round(output$P[round(output$P, digits = 3) != 0], digits = 3)
    output$print.p[round(output$P, digits = 3) == 0] <-
        format(output$P[round(output$P, digits = 3) == 0], digits = 2)
    output$sign <- sign(output$Coefficient)
    output$print.p <- sub("e", "*x*10^", output$print.p)
    if(use.ggplot){
      plot.bar(argv, output, prefix)
    }else{
      plot.bar.no.g(argv, output,prefix)
    }
}

plot.bar.no.g <- function(argv, output, prefix){
  png(paste(prefix, "_BARPLOT_", Sys.Date(), ".png", sep = ""),
      height=10, width=10, res=300, unit="in")
  layout(t(1:2), widths=c(8.8,1.2))
  par( cex.lab=1.5, cex.axis=1.25, font.lab=2, 
      oma=c(0,0.5,0,0),
      mar=c(4,6,0.5,0.5))
  xlab <- expression(italic(P) - value ~ threshold ~ (italic(P)[T]))
  ylab <- expression(paste("PRS model fit:  ", R ^ 2))
  if(argv$bar_col_p){
    col <- suppressWarnings(colorRampPalette(brewer.pal(12,argv$bar_palatte)))
    output <- output[order(output$Threshold),]
    output$color <-  col(nrow(output))
    b<- barplot(height=output$R2, 
                col=output$color, 
                border=NA, 
                ylim=c(0, max(output$R2)*1.25), 
                axes = F , ann=F)
    odd <- seq(0,nrow(output)+1,2)
    even <- seq(1,nrow(output),2)
    axis(side=1, at=b[odd], labels=output$Threshold[odd])
    axis(side=1, at=b[even], labels=output$Threshold[even])
    text( parse(text=paste(
      output$print.p)), 
      x = b+0.1, 
      y =  output$R2+ (max(output$R2)*1.05-max(output$R2)), 
      srt = 45)
  }else{
    col <- suppressWarnings(colorRampPalette(c(argv$bar_col_low, argv$bar_col_high)))
    output <- output[order(-log10(output$P)),]
    output$color <-  col(nrow(output))
    output <- output[order(output$Threshold),]
    b<- barplot(height=output$R2, 
                col=output$color, 
                border=NA, 
                ylim=c(0, max(output$R2)*1.25), 
                axes = F, ann=F)
    
    odd <- seq(0,nrow(output)+1,2)
    even <- seq(1,nrow(output),2)
    axis(side=1, at=b[odd], labels=output$Threshold[odd], lwd=2)
    axis(side=1, at=b[even], labels=output$Threshold[even],lwd=2)
    axis(side=1, at=c(0,b[1],2*b[length(b)]-b[length(b)-1]), labels=c("","",""), lwd=2, lwd.tick=0)
    text( parse(text=paste(
      output$print.p)), 
      x = b+0.1, 
      y =  output$R2+ (max(output$R2)*1.05-max(output$R2)), 
      srt = 45)
  }
  box(bty='L', lwd=2)
  axis(2,las=2, lwd=2)
  
  title(ylab=ylab, line=4, cex.lab=1.5, font=2 )
  title(xlab=xlab, line=2.5, cex.lab=1.5, font=2 )
  
  par(cex.lab=1.5, cex.axis=1.25, font.lab=2, 
      mar=c(20,0,20,4))
  if(argv$bar_col_p){
    output <- output[order(output$Threshold),]
    image(1, output$Threshold, t(output$Threshold), col=output$color, axes=FALSE, ann=F)
    axis(4,las=2,xaxs='r',yaxs='r', tck=0.2, col="white")
    title( bquote(atop(italic(P) - value , threshold), ),  
           line=2, cex=1.5, font=2, adj=0)
  }else{
    output <- output[order(-log10(output$P)),]
    image(1, -log10(output$P), t(seq_along(-log10(output$P))), col=output$color, axes=F,ann=F)
    axis(4,las=2,xaxs='r',yaxs='r', tck=0.2, col="white")
    title(bquote(atop(-log[10] ~ model, italic(P) - value), ), 
          line=2, cex=1.5, font=2, adj=0)
  }
  g<-dev.off()
}

plot.bar <- function(argv, output, prefix){
  ggfig.plot <- ggplot(data = output, aes(x = factor(Threshold), y = R2)) + geom_text(
    aes(label = paste(print.p)),
    vjust = -1.5,
    hjust = 0,
    angle = 45,
    cex = 4,
    parse = T
  )  +
    theme_sam + 
    scale_y_continuous(limits = c(0, max(output$R2) * 1.25)) +
    xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
    ylab(expression(paste("PRS model fit:  ", R ^ 2)))
  if (argv$bar_col_p) {
    ggfig.plot <-
      ggfig.plot + geom_bar(aes(fill = factor(Threshold)), stat = "identity") +
      scale_fill_brewer(palette = argv$bar_palatte,
                        name = expression(italic(P) - value ~ threshold))
  }else {
    ggfig.plot <-
      ggfig.plot + geom_bar(aes(fill = -log10(P)), stat = "identity") +
      scale_fill_gradient(
        low = argv$bar_col_low,
        high = argv$bar_col_high,
        name = bquote(atop(-log[10] ~ model, italic(P) - value), )
      )
  }
  
  ggsave(
    paste(prefix, "_BARPLOT_", Sys.Date(), ".png", sep = ""),
    ggfig.plot,
    height=10,width=10
  )
}


# Plot multi-phenotype plot -----------------------------------------------

multi_pheno_plot <- function(parameters, use.ggplot, use.data.table){
    writeLines("Plotting the Multi-Phenotype Plot")
    prs.summary <- NULL
    if(use.data.table){
        prs.summary <- fread(paste0(parameters$out, ".summary"), data.table=F)
    }else{
        prs.summary <- read.table(paste0(parameters$out, ".summary"), header=T)
    }
    multipheno <- subset(prs.summary, Set=="Base")
    multipheno$Phenotype <- sapply(multipheno$Phenotype, shorten_label)
    multipheno <- multipheno[order(multipheno$PRS.R2, decreasing=T), ]
    multipheno$Phenotype <- factor(multipheno$Phenotype, levels = multipheno$Phenotype)
    if(use.ggplot){
        b <-
            ggplot(multipheno[1:(min(parameters$multi_plot, nrow(multipheno))), ],
                   aes(
                       x = Phenotype,
                       y = PRS.R2,
                       fill = -log10(P)
                   )) +
            theme_sam +
            geom_bar(stat = "identity") +
            coord_flip() +
            ylab("Variance explained by PRS") +
            scale_fill_distiller(palette = "Spectral", name = bquote(atop(-log[10] ~ model, italic(P) - value), ))
        ggsave(paste(
            parameters$out,
            "_MULTIPHENO_BARPLOT_",
            Sys.Date(),
            ".png",
            sep = ""
        ))
    }else{
        png(paste(parameters$out, "_MULTIPHENO_BARPLOT_", Sys.Date(), ".png", sep = ""),
            height=10, width=10, res=300, unit="in")
        layout(t(1:2), widths=c(8.8,1.2))
        
        output <- multipheno[1:(min(parameters$multi_plot, nrow(multipheno))), ]
        max.label.length <- max(sapply(output$Phenotype, 
                                       max_length
        ))*0.75
        par( cex.lab=1.5, cex.axis=1.25, font.lab=2, 
             oma=c(0,0.5,0,0),
             mar=c(4,max.label.length,0.5,0.5))
        ylab <- "Variance explained by PRS"
        col <- suppressWarnings(colorRampPalette(brewer.pal(12,"RdYlBu")))
        
        output <- output[order(-log10(output$P)),]
        output$color <-  col(nrow(output))
        output <- output[order(output$PRS.R2),]
        b<- barplot(height=output$PRS.R2, 
                    col=output$color, 
                    border=NA, 
                    ann=F, horiz=TRUE,
                    ylab="",
                    xlab=ylab)
        axis(2,las=2, lwd=2, at=b, labels=output$Phenotype)
        box(bty="L",lwd=2)
        
        
        par(cex.lab=1.5, cex.axis=1.25, font.lab=2, 
            mar=c(20,0,20,4))
        output <- output[order(-log10(output$P)),]
        image(1, -log10(output$P), t(seq_along(-log10(output$P))), col=output$color, axes=F,ann=F)
        axis(4,las=2,xaxs='r',yaxs='r', tck=0.2, col="white")
        title(bquote(atop(-log[10] ~ model, italic(P) - value), ), 
              line=2, cex=1.5, font=2, adj=0)
        
        g<-dev.off()
    }
}
# Plot multi-set plot -----------------------------------------------------


multi_set_plot <- function(prefix, prs.summary, pheno.name, parameters, use.ggplot){
    writeLines("Plotting Multi-Set-Plot")
    if (nrow(prs.summary) < 1)
        stop((
            "Error: Cannot generate multi-plot as only one phenotype and the base set was observed!"
        ))
    overview <- subset(prs.summary, Phenotype==pheno.name)
    # process phenotype & pathway name to make it fit into the plot
    overview$Phenotype <- sapply(overview$Phenotype, shorten_label)
    overview$Set <- sapply(overview$Set, shorten_label)
    sets <- unique(overview$Set)
    multiset <- overview[order(overview$PRS.R2, decreasing=T), ]
    multiset$Set <- factor(multiset$Set, levels = multiset$Set)
    if(use.ggplot){
        b <-
            ggplot(multiset[1:(min(parameters$multi_plot, nrow(multiset))),], aes(
                x = Set,
                y = PRS.R2,
                fill = -log10(P)
            )) +
            theme_sam +
            geom_bar(stat = "identity") +
            coord_flip() +
            ylab("Variance explained by PRS") +
            scale_fill_distiller(palette = "PuOr", name = bquote(atop(-log[10] ~ model, italic(P) - value),)) +
            theme(axis.title.y = element_blank())
        ggsave(
            paste(prefix,
                  "_MULTISET_BARPLOT_",
                  Sys.Date(),
                  ".png",
                  sep = ""),
            b,
            height = 10,
            width = 10
        )
    } else{
        png(
            paste(
                prefix,
                "_MULTISET_BARPLOT_",
                Sys.Date(),
                ".png",
                sep = ""
            ),
            height = 10,
            width = 10,
            res = 300,
            unit = "in"
        )
        layout(t(1:2), widths = c(8.8, 1.2))
        output <-
            multiset[1:(min(parameters$multi_plot, nrow(multiset))),]
        max.label.length <- max(sapply(output$Set,
                                       max_length)) * 0.75
        par(
            cex.lab = 1.5,
            cex.axis = 1.25,
            font.lab = 2,
            oma = c(0, 0.5, 0, 0),
            mar = c(4, max.label.length, 0.5, 0.5)
        )
        ylab <- "Variance explained by PRS"
        col <-
            suppressWarnings(colorRampPalette(brewer.pal(12, "PuOr")))
        
        output <- output[order(-log10(output$P)), ]
        output$color <-  col(nrow(output))
        output <- output[order(output$PRS.R2, decreasing=T), ]
        b <- barplot(
            height = output$PRS.R2,
            col = output$color,
            border = NA,
            ann = F,
            horiz = TRUE,
            ylab = "",
            xlab = ylab
        )
        axis(
            2,
            las = 2,
            lwd = 2,
            at = b,
            labels = output$Set
        )
        box(bty = "L", lwd = 2)
        
        
        par(
            cex.lab = 1.5,
            cex.axis = 1.25,
            font.lab = 2,
            mar = c(20, 0, 20, 4)
        )
        output <- output[order(-log10(output$P)), ]
        image(
            1,
            -log10(output$P),
            t(seq_along(-log10(output$P))),
            col = output$color,
            axes = F,
            ann = F
        )
        axis(
            4,
            las = 2,
            xaxs = 'r',
            yaxs = 'r',
            tck = 0.2,
            col = "white"
        )
        title(
            bquote(atop(-log[10] ~ model, italic(P) - value),),
            line = 2,
            cex = 1.5,
            font = 2,
            adj = 0
        )
        
        g <- dev.off()
    }
}

# Sanity Check ------------------------------------------------------------


if (provided("no_regress", argv)) {
    quit("yes")
}
ignore_fid <- provided("ignore_fid", argv)

extract_matrix <- function(x, y) {
    z = which(x == y)
    
    if (length(z) == 0) {
        return(NA)
    } else{
        return(z)
    }
}
# CALL PLOTTING FUNCTION: Process the input names and call the actual plotting function
# we need to deduce the file names
# Now we actually require one single string for the input, separated by ,
# Get all the region information
# With this update, we only allow a single base file therefore we don't even need the
# information of base here

if (!provided("target", argv)) {
    stop("Target file name not found. You'll need to provide the target name for plotting! (even with --plot)")
}


phenos <- NULL

binary_target <- strsplit(argv$binary_target, split = ",")[[1]]
pheno.index <- 6
if (provided("pheno_col", argv)) {
    phenos <- strsplit(argv$pheno_col, split = ",")[[1]]
    if (!provided("pheno_file", argv)) {
        writeLines(
            strwrap(
                "WARNING: Cannot have multiple phenotypes if pheno_file is not provided. We will ignore the pheno_col option.",
                width = 80
            )
        )
    }else if (length(binary_target) != length(phenos)) {
        message <-
            "Number of binray target should match number of phenotype provided!"
        message <- paste(
            message,
            "There are ",
            length(binary_target),
            " binary target information and ",
            length(phenos),
            "phenotypes",
            sep = ""
        )
        stop(message)
    } else{
        header <- read.table(argv$pheno_file, nrows = 1, header = TRUE)
        # This will automatically filter out un-used phenos
        valid.pheno <- phenos %in% colnames(header)
        valid.file.index <- colnames(header) %in% phenos
        if (sum(valid.pheno) == 0) {
            stop("Error: None of the phenotype is identified in phenotype header!")
        }
        binary_target <- binary_target[valid.pheno]
        phenos <- phenos[valid.pheno]
        pheno.index <- c(1:ncol(header))[valid.file.index]
    }
} else if (provided("pheno_file", argv)) {
    pheno.index <- 3
    if (ignore_fid)
        pheno.index <- 2
} else{
    if (length(binary_target) != 1) {
        stop("Too many binary target information. We only have one phenotype")
    }
}


# Read in covariates ------------------------------------------------------
update_cov_header <- function(c) {
    res <- NULL
    for (i in c) {
        if (substr(i, 0, 1) == "@") {
            i <- substr(i, 2, nchar(i))
            temp <- strsplit(i, "\\[")[[1]]
            
            info <- NULL
            is_list <- NULL
            for (j in temp) {
                if (grepl("\\]", j)) {
                    tem <- strsplit(j, "\\]")[[1]]
                    
                    for (k in 1:length(tem)) {
                        info <- rbind(info, tem[k])
                        
                        is_list <- rbind(is_list, k == 1)
                        
                    }
                } else{
                    info <- rbind(info, j)
                    
                    is_list <- rbind(is_list, FALSE)
                    
                }
            }
            final <- NULL
            
            for (j in 1:nrow(info)) {
                if (is_list[j]) {
                    num <- NULL
                    ind <- strsplit(info[j], split = "\\.")[[1]]
                    
                    for (k in ind) {
                        if (grepl("-", k)) {
                            range <- strsplit(k , split = "-")[[1]]
                            
                            r <- range[1]:range[2]
                            
                            num <- c(num, r)
                            
                        } else{
                            num <- c(num, k)
                            
                        }
                    }
                    cur <- final
                    final <- NULL
                    for (n in num) {
                        final <- c(final, paste(cur, n, sep = ""))
                        
                    }
                } else{
                    final <- paste(final, info[j], sep = "")
                    
                }
            }
            res <- c(res, final)
        } else{
            res <- c(res, i)
            
        }
    }
    return(res)
}

covariance <- NULL
covariance.base <- NULL
if (provided("cov_file", argv)) {
    if (use.data.table) {
        covariance <- fread(argv$cov_file,
                            data.table = F,
                            header = T)
    } else {
        covariance <- read.table(argv$cov_file, header = T)
    }
    cov.header <- colnames(covariance)
    selected.cov <- cov.header[!cov.header%in%c("FID", "IID")]
    if(provided("cov_col", argv)){
        c <- strsplit(argv$cov_col, split = ",")[[1]]
        c <- update_cov_header(c)
        selected.cov <- cov.header[cov.header %in% c]
    }
    covariance.base <- covariance[, cov.header%in%c("FID", "IID",selected.cov)]
}

# we no longer have those complication
prefix <- argv$out

#regions <- read.table(paste(prefix, "region", sep = "."), header =T)
#num_region = nrow(regions)

# Process plot functions --------------------------------------------------

process_plot <-
    function(prefix,
             covariance,
             is_binary,
             pheno.file,
             parameters,
             pheno.index,
             use.data.table,
             use.ggplot,
             pheno.name) {
        sum.prefix <- prefix
        if(pheno.name!="-"){
            prefix <- paste(prefix, pheno.name, sep=".")
        }
        best <- NULL
        prs.summary <- NULL
        prsice.result <- NULL
        phenotype <- NULL
        if (use.data.table) {
            best <- fread(paste0(prefix, ".best"), data.table = F)
            prs.summary <-
                fread(paste0(sum.prefix, ".summary"), data.table = F)
            prsice.result <-
                fread(paste0(prefix, ".prsice"), data.table = F)
            phenotype <-
                fread(pheno.file, data.table = F, header = F)
        } else{
            best <- read.table(paste0(prefix, ".best"), header = T)
            prs.summary <-
                read.table(paste0(sum.prefix, ".summary"), header = T)
            prsice.result <-
                read.table(paste0(prefix, ".prsice"), header = T)
            # Allow header = false for fam or for phenotype files that does not contain phenotype name
            phenotype <- read.table(pheno.file, header = F)
        }
        best <- subset(best, In_Regression == "Yes")
        # We know the format of the best file, and it will always contain FID and IID
        
        base.prs <- best[,c(1,2,4)]
        if(provided("plot_set", parameters) & (provided("msigdb", parameters) | provided("bed", parameters) | provided("gtf", parameters)| provided("snp_set", parameters)| provided("snp_sets", parameters))){
            base.prs <- best[,colnames(best)%in%c("FID", "IID", parameters$plot_set)]
            colnames(base.prs)[3] <- "PRS"
        }
# Generate phenotype matrix -----------------------------------------------
        # extract the phenotype column
        # And only retain samples with phenotype and covariate information
        # They will be found in the best data.frame
        
        ignore_fid <- provided("ignore_fid", parameters)
        if (!ignore_fid) {
            phenotype <- phenotype[, c(1:2, pheno.index)]
            colnames(phenotype) <- c("FID", "IID", "Pheno")
            phenotype <-
                phenotype[phenotype$FID %in% best$FID &
                              phenotype$IID %in% best$IID, ]
        } else{
            phenotype <- phenotype[, c(1, pheno.index)]
            colnames(phenotype) <- c("IID", "Pheno")
            phenotype <- phenotype[phenotype$IID %in% best$IID, ]
        }
        phenotype$Pheno <- as.numeric(as.character(phenotype$Pheno))
        pheno <- phenotype
        use.residual <- F
        if(is_binary){
            if(max(pheno$Pheno)==2){
                pheno$Pheno <- pheno$Pheno-1
            }
        }
        if(!is.null(covariance)){
            # We will regress out the residual
            # Can direct merge as we have standardized the header
            temp.pheno <- merge(pheno, covariance)
            family <- gaussian
            if(is_binary){
                family <- binomial
            }
            residual <-
                rstandard(glm(Pheno ~ ., 
                              data = temp.pheno[, !colnames(temp.pheno) %in% c("FID", "IID")], 
                              family =family))
            pheno$Pheno <- residual
            use.residual <- T
        }
        
# Start calling functions -------------------------------------------------
        if (provided("quantile", parameters) && parameters$quantile > 0) {
            # Need to plot the quantile plot (Remember to remove the iid when performing the regression)
            if(!provided("quant_break", parameters)){
                quantile_plot(base.prs, pheno, prefix, parameters, is_binary, use.ggplot, use.residual)
            }else{
                uneven_quantile_plot(base.prs, pheno, prefix, parameters, is_binary, use.ggplot, use.residual)
            }
        }
        if(provided("msigdb", parameters) | provided("bed", parameters) | provided("gtf", parameters)
           | provided("snp_test", parameters) | provided("snp_tests", parameters)){
            if(length(strsplit(argv$bar_levels, split=",")[[1]])>1){
                bar_plot(prsice.result, prefix, parameters, use.ggplot) 
                if(!provided("fastscore", parameters)){
                    high_res_plot(prsice.result, prefix, parameters, use.ggplot)
                }
            }
        }else{
            bar_plot(prsice.result, prefix, parameters, use.ggplot)
            if(!provided("fastscore", parameters)){
                high_res_plot(prsice.result, prefix, parameters, use.ggplot)
            }
        }
        if(provided("multi_plot", parameters)){
            multi_set_plot(prefix, prs.summary, pheno.name, parameters, use.ggplot)
        }
    }


# Check if phenotype file is of sample format -----------------------------
is_sample_format <- function(file) {
    con = file(file, "r")
    first_line <- readLines(con, n = 1)
    second_line <- readLines(con, n = 1)
    close(con)
    if (length(first_line) == 0 | length(second_line) == 0) {
        # Unless there is only one sample? but that will be ridiculous for us to consider that
        stop("Error: Phenotype file should contain at least 2 line of input")
    }
    first <- strsplit(first_line, split = "\t")[[1]]
    if (length(first) == 1) {
        # Maybe seperated by space?
        first <- strsplit(first_line, split = " ")[[1]]
    }
    second <- strsplit(second_line, split = "\t")[[1]]
    if (length(first) == 1) {
        second <- strsplit(second_line, split = " ")[[1]]
    }
    if (length(first) != length(second) | length(first) < 3) {
        return(FALSE)
    }
    for (i in 1:3) {
        if (second[i] != 0) {
            return(FALSE)
        }
    }
    for (i in 4:length(second)) {
        if (second[i] != "D" &
            second[i] != "C" &
            second[i] != "P" & second[i] != "B") {
            return(FALSE)
        }
    }
    return(TRUE)
}
# Calling plot functions --------------------------------------------------
# Check target
pheno.file <- NULL

if (provided("pheno_file", argv)) {
    pheno.file <- argv$pheno_file
} else if(provided("target", argv)){
    # Check if external fam / sample file is provided
    target.info <- strsplit(argv$target, split = ",")[[1]]
    if (length(target.info) == 2) {
        pheno.file <- target.info[2]
        if (provided("type", argv)) {
            if (argv$type == "bgen") {
                # sample file should contain FID and IID by format requirement
                pheno.index <- 3
                if (ignore_fid &
                    !is_sample_format(pheno.file))
                    pheno.index <- 2
            }
        }
    } else{
        if (provided("type", argv)) {
            if (argv$type == "bgen") {
                stop("Error: You must provide either a phenotype or sample file for bgen input")
            } else if (argv$type == "bed") {
                pheno.file <- paste0(argv$target, ".fam")
            }
        } else{
            # Because default is always plink
            pheno.file <- paste0(argv$target, ".fam")
        }
    }
}else if(provided("target_list", argv)){
    # Assume no header
    target.list <- read.table(argv$target_list)
    target.prefix <- target.list[1]
    pheno.file <- paste0(target.prefix, ".fam")    
    if(provided("type", argv)){
        if(argv$type=="bgen"){
            stop("Error: You must provide a phenotype file for bgen list input")
        }   
    }
}

update_cov_factor <- function(parameters, pheno.file, pheno.index, cov.base){
    
    phenotype <- NULL
    if (use.data.table) {
        phenotype <-
            fread(pheno.file, data.table = F, header = F)
    } else{
        # Allow header = false for fam or for phenotype files that does not contain phenotype name
        phenotype <- read.table(pheno.file, header = F)
    }
    ignore_fid <- provided("ignore_fid", parameters)
    if (!ignore_fid) {
        phenotype <- phenotype[, c(1:2, pheno.index)]
        colnames(phenotype) <- c("FID", "IID", "Pheno")
        phenotype$Pheno <- suppressWarnings(as.numeric(as.character(phenotype$Pheno)))
        phenotype <-
            phenotype[!is.na(phenotype$Pheno), ]
    } else{
        phenotype <- phenotype[, c(1, pheno.index)]
        colnames(phenotype) <- c("IID", "Pheno")
        phenotype$Pheno <- suppressWarnings(as.numeric(as.character(phenotype$Pheno)))
        phenotype <-
            phenotype[!is.na(phenotype$Pheno), ]
    }
    # Now remove any missing sample from cov.base
    covariance <- NULL
    if(!is.null(cov.base)){
        if(!ignore_fid){
            covariance <- cov.base[cov.base$FID %in% phenotype$FID & 
                                       cov.base$IID %in% phenotype$IID, ]
        }else{
            covariance <- cov.base[cov.base$IID %in% phenotype$IID, ]
        }
        # Note: cov.base only contains the valid headers
        if(provided("cov_factor", parameters)){
            factor_cov <- parameters$cov_factor
            for(i in colnames(covariance)){
                if(i != "FID" & i != "IID"){
                    if(i %in%factor_cov){
                        covariance[,i] <- factor(covariance[,i], levels=unique(covariance[,i]))
                    }else{
                        covariance[,i] <- as.numeric(covariance[,i])
                    }
                }
            }
        }else{
            for(i in colnames(covariance)){
                if(i!="FID" & i!="IID"){
                    covariance[,i] <- as.numeric(covariance[,i])
                }
            }
        }
        if (ignore_fid) { 
            colnames(covariance) <- 
                c("IID", paste0("Cov", 1:(ncol(covariance) - 1))) 
        } else{ 
            colnames(covariance) <- 
                c("FID", "IID", paste0("Cov", 1:(ncol(covariance) - 2))) 
        } 
    }
    return(covariance)
}
# To account for the chromosome number
pheno.file <- gsub("#", "1", pheno.file)
if (!is.null(phenos) &
    length(phenos) > 1) {
    for (i in 1:length(phenos)) {
        # Update the covariance matrix accordingly
        covariance <- update_cov_factor(argv, pheno.file, pheno.index[i], covariance.base)
        process_plot(
            argv$out,
            covariance,
            binary_target[i],
            pheno.file,
            argv,
            pheno.index[i],
            use.data.table,
            use.ggplot,
            phenos[i]
        )
    }
    if(provided("multi_plot", argv)){
        multi_pheno_plot(argv, use.ggplot, use.data.table)
    }
} else if (!is.null(phenos)) {
    covariance <- update_cov_factor(argv, pheno.file, pheno.index[1], covariance.base)
    process_plot(
        argv$out,
        covariance,
        binary_target[1],
        pheno.file,
        argv,
        pheno.index[1],
        use.data.table,
        use.ggplot,
        "-"
    )
} else{
    covariance <- update_cov_factor(argv, pheno.file, pheno.index[1], covariance.base)
    process_plot(
        argv$out,
        covariance,
        binary_target[1],
        pheno.file,
        argv,
        pheno.index[1],
        use.data.table,
        use.ggplot,
        "-"
    )
}



