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
    --type                  File type of the target file. Support bed \n
                            (binary plink) and bgen format. Default: bed\n
\nDosage:\n
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
    --quant-extract | -e    File containing sample ID to be plot on a separated\n
                            quantile e.g. extra quantile containing only \n
                            schizophrenia samples. Must contain IID. Should\n
                            contain FID if --ignore-fid isn't set.\n
    --quant-ref             Reference quantile for quantile plot\n
    --scatter-r2            y-axis of the high resolution scatter plot should be R2\n

\nMisc:\n
    --all-score             Output PRS for ALL threshold. WARNING: This\n
                            will generate a huge file\n
    --exclude               File contains SNPs to be excluded from the\n
                            analysis\n
    --extract               File contains SNPs to be included in the \n
                            analysis\n
    --ignore-fid            Ignore FID for all input. When this is set,\n
                            first column of all file will be assume to\n
                            be IID instead of FID\n
    --logit-perm            When performing permutation, still use logistic\n
                            regression instead of linear regression. This\n
                            will substantially slow down PRSice\n
    --keep-ambig            Keep ambiguous SNPs. Only use this option\n
                            if you are certain that the base and target\n
                            has the same A1 and A2 alleles\n
    --out           | -o    Prefix for all file output\n
    --perm                  Number of permutation to perform. This swill\n
                            generate the empirical p-value. Recommend to\n
                            use value larger than 10,000\n
    --seed          | -s    Seed used for permutation. If not provided,\n
                            system time will be used as seed. When same\n
                            seed and same input is provided, same result\n
                            can be generated\n
    --print-snp             Print all SNPs used to construct the best PRS\n
    --thread        | -n    Number of thread use\n
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
      "RColorBrewer")
found <- FALSE
argv <- commandArgs(trailingOnly = TRUE)
dir_loc <- grep("--dir", argv)
if (length(dir_loc) != 0) {
    dir_loc <- dir_loc + 1
    found <- TRUE
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

UsePackage <- function(package, dir)
{
    if (!InstalledPackage(package))
    {
        dir.create(file.path(dir, "lib"), showWarnings = FALSE)
        .libPaths(c(.libPaths(), paste(dir, "/lib", sep = "")))
        if (!InstalledPackage(package)) {
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
    if (found)
    {
        if (!UsePackage(library, argv[dir_loc]))
        {
          if(library=="data.table"){
            use.data.table <- F
            writeLines("Cannot install data.table, will fall back and use read.table instead")
            writeLines("Note: It will be slower when reading large files")
          }else if(library=="ggplot2"){
            use.ggplot <- F
            writeLines("Cannot install ggplot2, will fall back and native plotting devices")
            writeLines("Note: The legends will be uglier")
          }else{
            stop("Error: ", library, " cannot be load nor install!")
          }
        }
    } else{
        if (!UsePackage(library, "."))
        {
            if(library=="data.table"){
              use.data.table <- F
              writeLines("Cannot install data.table, will fall back and use read.table instead")
              writeLines("Note: It will be slower when reading large files")
            }else if(library=="ggplot2"){
              use.ggplot <- F
              writeLines("Cannot install ggplot2, will fall back and native plotting devices")
              writeLines("Note: No legend will be displayed for bar-chart")
            }else{
              stop("Error: ", library, " cannot be load nor install!")
            }
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
  make_option(c("--type"), type = "character"),
  # Dosage
  make_option(c("--hard-thres"), type = "numeric"),
  make_option(c("--hard"), action = "store_true"),
  # Clumping
  make_option(c("--clump-kb"), type = "character", dest = "clump_kb"),
  make_option(c("--clump-r2"), type = "numeric", dest = "clump_r2"),
  make_option(c("--clump-p"), type = "numeric", dest = "clump_p"),
  make_option(c("-L", "--ld"), type = "character"),
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
  make_option(c("-n", "--thread"), type = "numeric"),
    #make_option(c("--no-x"), action = "store_true", dest = "no_x"),
    #make_option(c("--no-y"), action = "store_true", dest = "no_y"),
    #make_option(c("--no-xy"), action = "store_true", dest = "no_xy"),
    #make_option(c("--no-mt"), action = "store_true", dest = "no_mt"),
    #make_option(c("--num-auto"), type = "numeric", dest = "num_auto"),
    #R Specified options
    make_option(c("--plot"), action = "store_true"),
    make_option(c("--quantile", "-q"), type = "numeric"),
    make_option(c("--multi-plot"), type = "numeric", dest="multi_plot"),
    make_option(c("--plot-set"), type = "character", dest="plot_set", default="Base"),
    make_option(c("--quant-pheno"), action = "store_true", dest = "quant_pheno"),
    make_option(c("--quant-extract", "-e"), type = "character", dest = "quant_extract"),
    make_option("--quant-ref", type = "numeric", dest = "quant_ref"),
    make_option(
        "--scatter-r2",
        action = "store_true",
        default = F,
        dest = "scatter_r2"
    ),
    make_option(
        "--bar-col-p",
        action = "store_true",
        default = F,
        dest = "bar_col_p"
    ),
    make_option(
        "--bar-col-low",
        type = "character",
        default = "dodgerblue",
        dest = "bar_col_low"
    ),
    make_option(
        "--bar-col-high",
        type = "character",
        default = "firebrick",
        dest = "bar_col_high"
    ),
    make_option(
        "--bar-palatte",
        type = "character",
        default = "YlOrRd",
        dest = "bar_palatte"
    ),
    make_option("--prsice", type = "character"),
    make_option("--dir", type = "character")
)


capture <- commandArgs(trailingOnly = TRUE)
help <- (sum(c("--help", "-h") %in% capture) >= 1)
has_c <- (sum(c("--prsice") %in% capture) >= 1)
if (help && !has_c) {
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
    "scatter-r2",
    "bar-col-p",
    "bar-col-low",
    "bar-col-high",
    "bar-palatte",
    "prsice",
    "multi-plot",
    "plot-set",
    "dir"
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
# This is also one of the reason why window doesn't work. I don't know if we can handle the \
# WINDOW PEOPLE
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
        "beta",
        "no-full",
        "no-default",
        "ignore-fid",
        "index",
        "keep-ambig",
        "logit-perm",
        "no-clump",
        "no-regress",
        "no-x",
        "no-y",
        "no-xy",
        "no-mt",
        "fastscore",
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

# Read in the commands ----------------------------------------------------
logFile <- paste(argv$out,"log",sep=".")
con  <- file(logFile, open = "r")

c<- NULL
# Only need to know the information of binary target
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  line <- (trimws(oneLine))
  if(startsWith(line, "--")){
    commands <- strsplit(line, split=" ")[[1]]
    commands <- commands[commands!="" & commands!="\\"]
    if(length(commands) == 1){
      # Flag
    }else{
      # With input
      c <- gsub("--", "", commands[1])
      i <- commands[2]
      if(c=="binary-target"){
        argv$binary_target <- commands[2]
      }else if(c=="bar-levels"){
        argv$bar_levels <- commands[2]
      }
    }
  }
} 
close(con)
# Plottings ---------------------------------------------------------------

# Standard Theme for all plots

# Te selected theme to be used 
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
# PLOTTING: Here contains all the function for plotting
# quantile_plot: plotting the quantile plots
quantile_plot <-
    function(PRS, PRS.best, pheno, prefix, argv, binary, use.ggplot) {
        PRS.best <- PRS.best[,!colnames(PRS.best)%in%c("Has_Phenotype")]
        binary <- as.logical(binary)
        writeLines("Plotting the quantile plot")
        num_cov <- sum(!colnames(pheno)%in%c("FID", "IID", "Pheno"))
        extract = NULL
        if (provided("quant_extract", argv)) {
            if(use.data.table){
              extract <- fread(argv$quant_extract,
                            header = F,
                            data.table = F)
            }else{
              extract <- read.table(argv$quant_extract,
                              header = F)
            }
          
        }
        num_quant <- argv$quantile
        
        pheno.include <-
            NULL #Because we always name the phenotype as pheno, it will never be PRS
        
        if (provided("ignore_fid", argv) & sum(colnames(pheno)%in%c("FID"))!=1) {
            pheno.merge <- merge(PRS.best, pheno, by = "IID")
        } else{
            pheno.merge <- merge(PRS.best, pheno, by = c("FID", "IID"))
        }
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
        if(num_cov > 0){
          reg <- pheno.merge[,c("Pheno", paste("Cov", 1:num_cov))]
          family <- gaussian
          if(binary){
            family <- binomial 
          }
          residual <- rstandard(glm(Pheno~., family=family, data=reg))
          pheno.merge <- data.frame(Pheno=residual, PRS=pheno.merge$PRS)
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
        }
        if (!pheno.as.quant) {
          quants <- as.numeric(cut(
            pheno.merge$PRS,
            breaks = unique(quantile(
              pheno.merge$PRS, probs = seq(0, 1, 1 / num_quant)
            )),
            include.lowest = T
          ))
          
        } else{
          quants <- as.numeric(cut(
            pheno.merge$Pheno,
            breaks = unique(quantile(
              pheno.merge$Pheno, probs = seq(0, 1, 1 / num_quant)
            )),
            include.lowest = T
          ))
        }
        
        if (anyDuplicated(quantile(pheno.merge$PRS, probs = seq(0, 1, 1 / num_quant)))) {
            writeLines(paste(
                "Duplicate quantiles formed. Will use less quantiles: ",
                length(unique(quants)),
                sep = ""
            ))
            
        }
        num_quant = sum(!is.na(unique(quants)))
        if (!is.null(extract)) {
          extract_ID <- NULL
          best_ID <- NULL
          if(provided("ignore_fid",argv)){
            extract_ID <- extract$V1
            best_ID <- pheno.merge$IID
          }else{
            extract_ID <- paste(extract$V1, extract$V2, sep = "_")
            best_ID <- paste(pheno.merge$FID, pheno.merge$IID, sep = "_")
          }
          quants[best_ID %in% extract_ID] <-
            num_quant + 1 # We only matched based on the IID here
          num_quant <- num_quant + 1
        }
        quant.ref <- ceiling(argv$quantile / 2)
        if (provided("quant_ref", argv)) {
          quant.ref <- argv$quant_ref
          if (quant.ref > argv$quantile) {
            quant.ref <- ceiling(argv$quantile / 2)
            writeLines(
              paste(
                "WARNING: reference quantile",
                quant.ref,
                "is greater than number of quantiles",
                argv$quantile,
                "\n Using middle quantile by default"
              )
            )
          }
        }
        
        quants <-
          factor(quants, levels = c(quant.ref, seq(min(quants), max(quants), 1)[-quant.ref]))
        pheno.merge$quantile <- quants
        if (!pheno.as.quant) {
            family <- gaussian
            if (binary) {
                if( num_cov == 0 ){
                    family <- binomial
                }
            }
            reg <- summary(glm(Pheno ~ quantile, family, data = pheno.merge))
            coef.quantiles <- (reg$coefficients[1:num_quant, 1])
            ci <- (1.96 * reg$coefficients[1:num_quant, 2])
            
            ci.quantiles.u <-
                coef.quantiles + ci
            ci.quantiles.l <-
                coef.quantiles - ci
            if(binary){
                if(num_cov==0){
                    ci.quantiles.u <- exp(ci.quantiles.u)
                    ci.quantiles.l <- exp(ci.quantiles.l)
                    coef.quantiles <- exp(coef.quantiles)
                }
            }
            
            coef.quantiles[1] <- ifelse(binary & (num_cov==0),1,0)
            ci.quantiles.u[1] <- ifelse(binary & (num_cov==0),1,0)
            ci.quantiles.l[1] <- ifelse(binary & (num_cov==0),1,0)
            quantiles.for.table <-
                c(quant.ref, seq(1, num_quant, 1)[-quant.ref])
            quantiles.df <-
                data.frame(
                    Coef = coef.quantiles,
                    CI.U = ci.quantiles.u,
                    CI.L = ci.quantiles.l,
                    DEC = quantiles.for.table
                )
            quantiles.df$Group = 0
            if (!is.null(extract)) {
                # Because the last quantile is set to be cases
                quantiles.df$Group[max(quantiles.df$DEC)] = 1
            }
            quantiles.df$Group <-
                factor(quantiles.df$Group, levels = c(0, 1))
            quantiles.df <- quantiles.df[order(quantiles.df$DEC), ]
            
            if(use.ggplot){
              plot.quant(quantiles.df, num_quant, binary, extract, prefix, num_cov)
            }else{
              plot.quant.no.g(quantiles.df, num_quant, binary, extract, prefix, num_cov)
            }
        }else{
            # TODO: Maybe also change this to regression? Though might be problematic if we have binary pheno without cov
            pheno.sum <- data.frame(mean=numeric(num_quant), quantile=1:num_quant, UCI=numeric(num_quant), LCI=numeric(num_quant))
            for(i in 1:num_quant){
                cur.prs <- pheno.merge$PRS[as.numeric(as.character(pheno.merge$quantile))%in%i]
                pheno.sum$mean[i] <-mean(cur.prs,na.rm=T)
                pheno.sum$UCI[i] <- pheno.sum$mean[i]+sd(cur.prs,na.rm=T)
                pheno.sum$LCI[i] <- pheno.sum$mean[i]-sd(cur.prs,na.rm=T)
            }
            pheno.sum$Group = 0
            if (!is.null(extract)) {
                pheno.sum$Group[num_quant] = 1
            }
            pheno.sum$Group <-
                factor(pheno.sum$Group, levels = c(0, 1))
            if(use.ggplot){
              plot.pheno.quant(pheno.sum, num_cov, num_quant, extract, prefix)
            }else{
              plot.pheno.quant.no.g(pheno.sum, num_cov, num_quant, extract, prefix)
            }
        }
        
    }

plot.pheno.quant.no.g <- function(pheno.sum, num_cov, num_quant, extract, prefix){
  png(paste(prefix, "_QUANTILES_PHENO_PLOT_", Sys.Date(), ".png", sep = ""),
      height=10, width=10, res=300, unit="in")
  par(pty="s", cex.lab=1.5, cex.axis=1.25, font.lab=2, mai=c(0.5,1.25,0.1,0.1))
  pheno.sum$color <- "#D55E00"
  xlab <- NULL
  if(num_cov>0){
    xlab <-"Quantiles for Residualized Phenotype"
  }else{
    xlab <-"Quantiles for Phenotype"
  }
  
  
  if(!is.null(extract)){
    pheno.sum$color <- "#D55E00"
    pheno.sum$color[quantiles.df$Group==1] <- "#0072B2"
  }
  ylab <- "Mean PRS given phenotype in quantiles"
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

plot.pheno.quant <- function(pheno.sum, num_cov, num_quant, extract, prefix){
  quantiles.plot <-
    ggplot(pheno.sum, aes(
      x = quantile,
      y = mean,
      ymin = LCI,
      ymax = UCI
    ))+ 
    theme_sam+
    scale_x_continuous(breaks = seq(0, num_quant, 1))+
    ylab("Mean PRS given phenotype in quantiles")
  if(num_cov>0){
    quantiles.plot <- quantiles.plot+
      xlab("Quantiles for Residualized Phenotype")
  }else{
    quantiles.plot <- quantiles.plot+
      xlab("Quantiles for Phenotype")
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

plot.quant <- function(quantiles.df, num_quant, binary, extract, prefix, num_cov){
    quantiles.plot <-
        ggplot(quantiles.df, aes(
            x = DEC,
            y = Coef,
            ymin = CI.L,
            ymax = CI.U
        )) + 
        theme_sam+
        scale_x_continuous(breaks = seq(0, num_quant, 1))+
        xlab("Quantiles for Polygenic Score")
    if (binary){
        if(num_cov==0) {
            quantiles.plot <-
                quantiles.plot + ylab("Odds Ratio for Score on Phenotype")
        }
    } else if(num_cov!=0){
        quantiles.plot <- quantiles.plot + ylab("Change in residualized\nPhenotype given score in quantiles")
    }else{
        quantiles.plot <- quantiles.plot +
            ylab("Change in Phenotype \ngiven score in quantiles")
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

plot.quant.no.g <- function(quantiles.df, num_quant, binary, extract, prefix, num_cov){
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
      if(num_cov==0) {
          quantiles.plot <-
              quantiles.plot + ylab("Odds Ratio for Score on Phenotype")
      }
  } else if(num_cov!=0){
      quantiles.plot <- quantiles.plot + ylab("Change in residualized\nPhenotype given score in quantiles")
  }else{
      quantiles.plot <- quantiles.plot +
          ylab("Change in Phenotype \ngiven score in quantiles")
  }
  xlab <- "Quantiles for Polygenic Score"
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

# run_plot: The function used for calling different plotting functions
run_plot <- function(prefix, argv, pheno_matrix, binary) {
    writeLines("")
    PRS <- NULL
    PRS.best <- NULL
    if(use.data.table){
      PRS <-
        fread(paste(prefix, ".prsice", sep = ""),
              data.table = F)
      PRS.best <-
        fread(paste(prefix, ".best", sep = ""),
              data.table = F)
    }else{
      
      PRS <-
        read.table(paste(prefix, ".prsice", sep = ""),
              header = T)
      PRS.best <-
        read.table(paste(prefix, ".best", sep = ""),
              header = T)
    }
    
    PRS.best <- subset(PRS.best, PRS.best$Has_Phenotype == "Yes")
    # colnames(PRS.best)[3] <- "PRS"
    # start from here, we need to organize all the file accordingly so that the individual actually match up with each other
    # Good thing is, only quantile plot really needs the cov and phenotype information
    if (provided("quantile", argv) && argv$quantile > 0) {
        # Need to plot the quantile plot (Remember to remove the iid when performing the regression)
        quantile_plot(PRS, PRS.best, pheno_matrix, prefix, argv, binary, use.ggplot)
      
    }
    # Now perform the barplotting
    if (!provided("fastscore", argv) || !argv$fastscore) {
        high_res_plot(PRS, prefix, argv, use.ggplot)
    }
    bar_plot(PRS, prefix, argv, use.ggplot)
}




# Process file names for plotting------------------------------------------------------
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
    stop("Target file name not found. You'll need to provide the target name for plotting!")
}

if (!provided("pheno_file", argv) &&
    !provided("binary_target", argv)) {
    argv$binary_target = "T"
    
} else{
    if (!provided("pheno_col", argv) && provided("binary_target", argv)
        && length(argv$binary_target) == 1) {
        #This is ok
        test <- "ok"
    } else if (!provided("pheno_col", argv) &&
               !provided("binary_target", argv)) {
        argv$binary_target = "T"
        
    } else if (provided("pheno_col", argv) &&
               !provided("binary_target", argv)
               && length(argv$pheno_col) <= 1) {
        argv$binary_target = "T"
    }
    else if (provided("pheno_col", argv) &&
             provided("binary_target", argv)
             && length(argv$pheno_col) != length(argv$binary_target)) {
        stop(
            "ERROR: Number of target phenotypes doesn't match information of binary target!
            You must indicate whether the phenotype is binary using --binary-target"
        )
    }
    }

phenos = NULL
binary_target = strsplit(argv$binary_target, split = ",")[[1]]
if (provided("pheno_col", argv)) {
    phenos = strsplit(argv$pheno_col, split = ",")[[1]]
    if (!provided("pheno_file", argv)) {
        writeLines(
            strwrap(
                "WARNING: Cannot have multiple phenotypes if pheno_file is not provided. We will ignore the pheno_col option.",
                width = 80
            )
        )
    } else{
        header <- read.table(argv$pheno_file, nrows = 1, header = TRUE)
        # This will automatically filter out un-used phenos
        if (length(binary_target) != length(phenos)) {
            message <-
                "Number of binray target should match number of phenotype provided!"
            message = paste(
                message,
                "There are ",
                length(binary_target),
                " binary target information and ",
                length(phenos),
                "phenotypes",
                sep = ""
            )
            stop(message)
        }
        binary_target = binary_target[phenos %in% colnames(header)]
        phenos = phenos[phenos %in% colnames(header)]
    }
} else{
    if (length(binary_target) != 1) {
        stop("Too many binary target information. We only have one phenotype")
    }
}



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

# Now we have the correct header of phenos, bases and binary_target information
# Need to get the covariates
covariance = NULL
if (provided("cov_file", argv)) {
    covariance <- NULL
    if(use.data.table){
      covariance <- fread(argv$cov_file, data.table = F, header = T)
    }else {
      covariance <- read.table(argv$cov_file, header = T)
    }
      
    if (provided("cov_col", argv)) {
        c = strsplit(argv$cov_col, split = ",")[[1]]
        c <- update_cov_header(c)
        selected <- colnames(covariance) %in% c
        if (!ignore_fid) {
            selected[2] <-
                TRUE #When ignore_fid isn't provided, then we need to also include the FID information
        }
        selected[1] <-
            TRUE # we always want the IID information, which should be the first column
        covariance <- covariance[, selected]
    }
    if (ignore_fid) {
        colnames(covariance) <- c("IID", paste("Cov", 1:(ncol(covariance) - 1)))
    } else{
        colnames(covariance) <-
            c("FID", "IID", paste("Cov", 1:(ncol(covariance) - 2)))
    }
}

# we no longer have those complication
prefix <- argv$out


regions <- read.table(paste(prefix, "region", sep = "."), header =
                          T)
num_region = nrow(regions)


region = argv$plot_set
# Do this for each phenotype

#fam = fread(paste(argv$target, ".fam", sep = ""), data.table = F, header = F )
#colnames(fam)[1:2] = c("FID", "IID")
#fam_id <- paste(fam$FID, fam$IID, sep="_")
#match_cov = NULL
#if (!is.null(covariance)) {
#  if(ignore_fid){
#    match_cov <- covariance[covariance$IID %in% fam$IID,]
#    match_cov <- match_cov[match(fam$IID, match_cov$IID),] # match the ordering
#    match_cov <- match_cov[ !is.na(apply(match_cov[,-1],1,sum)),]
#  }else{
#    cov_ID <- paste(covariance$FID, covariance$IID, sep="_")
#    match_cov = covariance[cov_ID %in% fam_id,]
#    cov_ID <- paste(match_cov$FID, match_cov$IID, sep="_") #updated cov_ID
#    match_cov = match_cov[match(fam_id, cov_ID),] # match the ordering
#    match_cov <- match_cov[ !is.na(apply(match_cov[,-c(1:2)],1,sum)),] #remove all NA
#  }
#}
#Now match_cov contain all samples with valid covariates


if (!is.null(phenos)) {
    pheno_file <- NULL
    if(use.data.table){
      pheno_file <- fread(argv$pheno_file,
                       data.table = F)
    }else{
      pheno_file <- read.table(argv$pheno_file,
                         header = T)
    }
    if (ignore_fid) {
        colnames(pheno_file)[1] <- "IID"
    } else{
        colnames(pheno_file)[1:2] <- c("FID", "IID")
    }
    id = 1
    phenos.index <-
        unlist(apply(
            as.matrix(phenos),
            1,
            extract_matrix ,
            colnames(pheno_file)
        ))
    
    for (p in 1:length(phenos.index)) {
        if (!is.na(phenos.index[p])) {
            cur_prefix <- prefix
            if (length(phenos.index) != 1) {
                cur_prefix <- paste(cur_prefix, phenos[p], sep = ".")
            }
            if (num_region != 1) {
                cur_prefix = paste(cur_prefix, region, sep = ".")
            }
            # Get the best score
            best<-NULL
            if(use.data.table){
            best <-
                fread(
                    paste(cur_prefix, "best", sep = "."),
                    data.table = F,
                    header = T
                )
            }else{
              
              best <-
                read(
                  paste(cur_prefix, "best", sep = "."),
                  header = T
                )
            }
            # Give run_plot a ready to use matrix
            cur_pheno <- NULL
            if (ignore_fid) {
                cur_pheno <-
                    data.frame(IID = pheno_file[, 1], Pheno = pheno_file[phenos.index[p]])
                colnames(cur_pheno)[2] <- "Pheno"
            } else{
                cur_pheno <-
                    data.frame(FID = pheno_file[, 1],
                               IID = pheno_file[, 2],
                               Pheno = pheno_file[phenos.index[p]])
                colnames(cur_pheno)[3] <- "Pheno"
            }
            if (binary_target[id]) {
                # Update the cur_pheno
                cur_pheno$Pheno <-
                    suppressWarnings(as.numeric(as.character(cur_pheno$Pheno)))
                cur_pheno$Pheno[cur_pheno$Pheno > 2] <- NA
                cur_pheno$Pheno[cur_pheno$Pheno < 0] <- NA
                if (max(cur_pheno$Pheno, na.rm = T) == 2 &&
                    min(cur_pheno$Pheno, na.rm = T) == 0) {
                    stop("Invalid case control formating. Either use 0/1 or 1/2 coding")
                } else if (max(cur_pheno$Pheno, na.rm = T) == 2) {
                    cur_pheno$Pheno = cur_pheno$Pheno - 1
                }
            }
            if (ignore_fid) {
                cur_pheno <- cur_pheno[cur_pheno$IID %in% best$IID, ]
            } else{
                cur_id <- paste(cur_pheno$FID, cur_pheno$IID, sep = "_")
                best_id <- paste(best$FID, best$IID, sep = "_")
                cur_pheno <- cur_pheno[cur_id %in% best_id, ]
            }
            cur_pheno.cov <- cur_pheno
            if (!is.null(covariance)) {
                if (ignore_fid) {
                    cur_pheno.cov = merge(cur_pheno, covariance, by = "IID")
                    
                } else{
                    cur_pheno.cov = merge(cur_pheno, covariance, by = c("FID", "IID"))
                    
                }
            }
            # This is a potential bug: What if there is duplicated IID?
            cur_pheno.cov = cur_pheno.cov[match(cur_pheno$IID, cur_pheno.cov$IID), ]
            run_plot(cur_prefix, argv, cur_pheno.cov, binary_target[id])
        } else{
            writeLines(paste(
                phenos[p],
                "not found in the phenotype file. It will be ignored"
            ))
        }
        id = id + 1
    }
} else{
    # No phenotype headers
    cur_prefix = paste(prefix, region, sep = ".")
    if (num_region == 1) {
        cur_prefix = paste(prefix, sep = ".")
    }
    pheno = NULL
    if (provided("pheno_file", argv)) {
        pheno <- NULL
        if(use.data.table){
          pheno <- fread(paste(argv$pheno_file),
                      data.table = F,
                      header = F)
        }else{
          pheno <- read.table(paste(argv$pheno_file),
                        header = F)
        }
        if (ignore_fid) {
            colnames(pheno)[1] <- "IID"
        } else{
            colnames(pheno)[1:3] <- c("FID", "IID", "V2") #Otherwise this is V3
        }
        
    }
    # Give run_plot a ready to use matrix
    best <- NULL
    if(use.data.table){
      best <-
        fread(paste(cur_prefix, "best", sep = "."),
              data.table = F,
              header = T)
    }else{
      best <-
        read.table(paste(cur_prefix, "best", sep = "."),
              header = T)
    }
    fam.clean = NULL
    if (!is.null(pheno))
    {
        if (ignore_fid) {
            fam.clean <- data.frame(IID = pheno$IID, Pheno = pheno$V2)
        } else{
            fam.clean <-
                data.frame(
                    FID = pheno$FID,
                    IID = pheno$IID,
                    Pheno = pheno$V2
                )
        }
    } else{
        fam_name <- argv$target
        fam_name <- gsub("#", "1", fam_name)
        fam <- NULL
        if(use.data.table){
          fam <-
            fread(paste(fam_name, "fam", sep = "."),
                  data.table = F,
                  header = F)
        }else{
          fam <-
            read.table(paste(fam_name, "fam", sep = "."),
                  header = F)
        }
        fam.clean <-
            data.frame(FID = fam$V1,
                       IID = fam$V2,
                       Pheno = fam$V6)
    }
    if (binary_target[1]) {
        # Update the cur_pheno
        fam.clean$Pheno <-
            suppressWarnings(as.numeric(as.character(fam.clean$Pheno)))
        fam.clean$Pheno[fam.clean$Pheno > 2] <- NA
        fam.clean$Pheno[fam.clean$Pheno < 0] <- NA
        if (max(fam.clean$Pheno, na.rm = T) == 2 &&
            min(fam.clean$Pheno, na.rm = T) == 0) {
            stop("Invalid case control formating. Either use 0/1 or 1/2 coding")
        } else if (max(fam.clean$Pheno, na.rm = T) == 2) {
            fam.clean$Pheno = fam.clean$Pheno - 1
        }
    }
    if (ignore_fid) {
        fam.clean <- fam.clean[fam.clean$IID %in% best$IID, ]
    } else{
        fam_id <- paste(fam.clean$FID, fam.clean$IID, "_")
        best_id <- paste(best$FID, best$IID, "_")
        fam.clean <- fam.clean[fam_id %in% best_id,]
    }
    cur_pheno.clean <- fam.clean
    fam.final <- cur_pheno.clean
    if (!is.null(covariance)) {
        if (ignore_fid) {
            fam.final <- merge(cur_pheno.clean, covariance, by = "IID")
        } else{
            fam.final <- merge(cur_pheno.clean, covariance, by = c("FID", "IID"))
        }
    }
    #Again, this can be error prone
    fam.final <-
        fam.final[match(cur_pheno.clean$IID, fam.final$IID), ]
    run_plot(cur_prefix, argv, fam.final, binary_target[1])
}


# Now check if the overview file is present
if (provided("multi_plot", argv)) {
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
    overview.name <- paste(argv$out, ".summary", sep = "")
    if (file.exists(overview.name)) {
        
        writeLines("Plotting Multi-Plot")
        overview <- read.table(overview.name, header = T)
        if (nrow(overview) < 1)
            stop((
                "Error: Cannot generate multi-plot as only one phenotype and the base set was observed!"
            )
            )
        overview$Phenotype <- sapply(overview$Phenotype, shorten_label)
        overview$Set <- sapply(overview$Set, shorten_label)
        phenos <- unique(overview$Phenotype)
        sets <- unique(overview$Set)
        if (length(phenos) != 1) {
            multipheno <- subset(overview, Set == "Base")
            multipheno <- multipheno[order(multipheno$PRS.R2), ]
            multipheno$Phenotype <-
                factor(multipheno$Phenotype, levels = multipheno$Phenotype)
            if(use.ggplot){
              b <-
                  ggplot(multipheno[1:(min(argv$multi_plot, nrow(multipheno))), ],
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
                  argv$out,
                  "_MULTIPHENO_BARPLOT_",
                  Sys.Date(),
                  ".png",
                  sep = ""
              ))
            }else{
                png(paste(argv$out, "_MULTIPHENO_BARPLOT_", Sys.Date(), ".png", sep = ""),
                    height=10, width=10, res=300, unit="in")
                layout(t(1:2), widths=c(8.8,1.2))
                
                output <- multipheno[1:(min(argv$multi_plot, nrow(multipheno))), ]
                max.label.length <- max(sapply(output$Phenotype, 
                                               function(x){
                                                   info <- strsplit(as.character(x), split="\n")[[1]];
                                                   max(sapply(info, nchar))
                                               }
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
            for (p in phenos) {
                multiset <- subset(overview, Phenotype == p)
                multiset <- multiset[order(multiset$PRS.R2), ]
                multiset$Set <-
                    factor(multiset$Set, levels = multiset$Set)
                if(use.ggplot){
                  b <-
                      ggplot(multiset[1:(min(argv$multi_plot, nrow(multiset))), ], aes(
                          x = Set,
                          y = PRS.R2,
                          fill = -log10(P)
                      )) +
                      theme_sam +
                      geom_bar(stat = "identity") +
                      coord_flip() +
                      ylab("Variance explained by PRS") +
                      scale_fill_distiller(palette = "PuOr", name = bquote(atop(-log[10] ~ model, italic(P) - value), )) +
                      theme(axis.title.y = element_blank())
                  ggsave(
                      paste(
                          argv$out,
                          "_",
                          p,
                          "_MULTISET_BARPLOT_",
                          Sys.Date(),
                          ".png",
                          sep = ""
                      ),
                      b,
                      height = 10,
                      width = 10
                  )
                }else{
                    png(paste(argv$out,
                              "_",
                              p,
                              "_MULTISET_BARPLOT_",
                              Sys.Date(),
                              ".png",
                              sep = ""
                    ),
                        height=10, width=10, res=300, unit="in")
                    layout(t(1:2), widths=c(8.8,1.2))
                    output <- multiset[1:(min(argv$multi_plot, nrow(multiset))), ]
                    max.label.length <- max(sapply(output$Set, 
                                                   function(x){
                                                       info <- strsplit(as.character(x), split="\n")[[1]];
                                                       max(sapply(info, nchar))
                                                    }
                                                   ))*0.75
                    par( cex.lab=1.5, cex.axis=1.25, font.lab=2, 
                         oma=c(0,0.5,0,0),
                         mar=c(4,max.label.length,0.5,0.5))
                    ylab <- "Variance explained by PRS"
                    col <- suppressWarnings(colorRampPalette(brewer.pal(12,"PuOr")))
                    
                    output <- output[order(-log10(output$P)),]
                    output$color <-  col(nrow(output))
                    output <- output[order(output$PRS.R2),]
                    b<- barplot(height=output$PRS.R2, 
                                col=output$color, 
                                border=NA, 
                                ann=F, horiz=TRUE,
                                ylab="",
                                xlab=ylab)
                    axis(2,las=2, lwd=2, at=b, labels=output$Set)
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
        } else{
            # Only plot one set plot. If phenotype == "-", replace it with pheno
            multiset <- overview
            multiset <- multiset[order(multiset$PRS.R2), ]
            multiset$Set <-
                factor(multiset$Set, levels = multiset$Set)
            if(use.ggplot){
              b <-
                  ggplot(multiset[1:(min(argv$multi_plot, nrow(multiset))), ], aes(
                      x = Set,
                      y = PRS.R2,
                      fill = -log10(P)
                  )) +
                  theme_sam +
                  geom_bar(stat = "identity") +
                  coord_flip() +
                  ylab("Variance explained by PRS") +
                  scale_fill_distiller(palette = "PuOr", name = bquote(atop(-log[10] ~ model, italic(P) - value), )) +
                  theme(axis.title.y = element_blank())
              ggsave(
                  paste(
                      argv$out,
                      "_MULTISET_BARPLOT_",
                      Sys.Date(),
                      ".png",
                      sep = ""
                  ),
                  b,
                  height = 10,
                  width = 10
              )
            }else{
                png(paste(argv$out,
                          "_MULTISET_BARPLOT_",
                          Sys.Date(),
                          ".png",
                          sep = ""
                    ),
                height=10, width=10, res=300, unit="in")
                layout(t(1:2), widths=c(8.8,1.2))
                output <- multiset[1:(min(argv$multi_plot, nrow(multiset))), ]
                max.label.length <- max(sapply(output$Set, 
                                               function(x){
                                                   info <- strsplit(as.character(x), split="\n")[[1]];
                                                   max(sapply(info, nchar))
                                               }
                ))*0.75
                par( cex.lab=1.5, cex.axis=1.25, font.lab=2, 
                     oma=c(0,0.5,0,0),
                     mar=c(4,max.label.length,0.5,0.5))
                ylab <- "Variance explained by PRS"
                col <- suppressWarnings(colorRampPalette(brewer.pal(12,"PuOr")))
                
                output <- output[order(-log10(output$P)),]
                output$color <-  col(nrow(output))
                output <- output[order(output$PRS.R2),]
                b<- barplot(height=output$PRS.R2, 
                            col=output$color, 
                            border=NA, 
                            ann=F, horiz=TRUE,
                            ylab="",
                            xlab=ylab)
                axis(2,las=2, lwd=2, at=b, labels=output$Set)
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
        
        
    }
}
