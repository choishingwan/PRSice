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
help_message <-"usage: Rscript PRSice.R [options] <-b base_file> <-t target_file> <--prsice prsice_location>\n
\nRequired:\n
    --prsice                Location of the PRSice binary\n
    --dir                   Location to install ggplot. Only require if ggplot\n
                            is not installed\n
\nBase File:\n
    --base          | -b    Base association file\n
    --beta                  Whether the test statistic is in the form of \n
                            BETA or OR. If set, test statistic is assume\n
                            to be in the form of BETA.\n
    --A1                    Column header containing the reference allele\n
                            Default: A1\n
    --A2                    Column header containing the alternative allele\n
                            Default: A2\n
    --bp                    Column header containing the SNP coordinate\n
                            Default: BP\n
    --chr                   Column header containing the chromosome\n
                            Default: CHR\n
    --index                 If set, assume the INDEX instead of NAME of\n
                            the corresponding columns are provided. Index\n
                            should be 0-based (start counting from 0)\n
    --pvalue        | -p    Column header containing the p-value\n
                            Default: P\n
    --se                    Column header containing the standard error\n
                            Default: SE\n
    --snp                   Column header containing the SNP ID\n
                            Default: SNP\n
    --stat                  Column header containing the summary statistic\n
                            If --beta is set, default as BETA. Otherwise,\n
                            try and search for OR or BETA from the header\n
                            of the base file\n
\nClumping:\n
    --clump-kb              The distance for clumping in kb\n
                            Default: 250\n
    --clump-r2              The R2 threshold for clumping\n
                            Default: 0.1\n
    --clump-p               The p-value threshold use for clumping.\n
                            Default: 1\n
    --ld            | -L    LD reference file. Use for LD calculation. If not\n
                            provided, will use the post-filtered target genotype\n
                            for LD calculation. Support multiple chromosome input\n
                            Please see --target for more information\n
    --ld-keep               File containing the sample(s) to be extracted from\n
                            the LD reference file. First column should be FID and\n
                            the second column should be IID. If --ignore-fid is\n
                            set, first column should be IID\n
                            Mutually exclusive from --ld-remove\n
    --ld-remove             File containing the sample(s) to be removed from\n
                            the LD reference file. First column should be FID and\n
                            the second column should be IID. If --ignore-fid is\n
                            set, first column should be IID\n
                            Mutually exclusive from --ld-keep\n
    --ld-type               File type of the LD file. Support bed (binary plink)\n
                            and bgen format. Default: bed\n
    --no-clump              Avoid performing clumping\n
    --proxy                 Proxy threshold for index SNP to be considered\n
                            as part of the region represented by the clumped\n
                            SNP(s). e.g. --proxy 0.8 means the index SNP will\n
                            represent region of any clumped SNP(s) that has a\n
                            R2>=0.8 even if the index SNP does not physically\n
                            locate within the region\n
\nCovariate:\n
    --cov-file      | -C    Covariate file. First column should be FID and \n
                            the second column should be IID. If --ignore-fid\n
                            is set, first column should be IID\n
    --cov-header    | -c    Header of covariates. If not provided, will use\n
                            all variables in the covariate file\n
\nDosage:\n
    --hard-thres            Hard threshold for dosage data. Any call less than\n
                            this will be treated as missing. Note that if dosage\n
                            data, is used as a LD reference, it will always be\n
                            hard coded to calculate the LD\n
                            Default: 0.9\n
    --hard                  Use hard coding instead of dosage for PRS construction.\n
                            Default is to use dosage instead of hard coding\n
\nPlotting:\n
    --bar-col-high          Colour of the most predicting threshold\n
                            Default: firebrick\n
    --bar-col-lower         Colour of the poorest predicting threshold\n
                            Default: dodgerblue\n
    --bar-col-p             Change the colour of bar to p-value threshold\n
                            instead of the association with phenotype\n
    --bar-palatte           Colour palatte to be used for bar plotting when\n
                            --bar_col_p is set. Default: YlOrRd\n
    --plot                  When set, will only perform plotting.\n
    --quantile      | -q    Number of quantiles to plot. No quantile plot\n
                            will be generated when this is not provided.\n
    --quant-extract | -e    File containing sample ID to be plot on a separated\n
                            quantile e.g. extra quantile containing only \n
                            schizophrenia samples. Must contain IID. Should\n
                            contain FID if --ignore-fid isn't set.\n
    --quant-ref             Reference quantile for quantile plot\n
    --scatter-r2            y-axis of the high resolution scatter plot should be R2\n
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
\nPRSice:\n
    --bar-levels            Level of barchart to be plotted. When --fastscore\n
                            is set, PRSice will only calculate the PRS for \n
                            threshold within the bar level. Levels should be\n
                            comma separated without space\n
    --fastscore             Only calculate threshold stated in --bar-levels\n
    --full                  Include the full model in the analysis\n
    --interval      | -i    The step size of the threshold. Default: 0.00005\n
    --lower         | -l    The starting p-value threshold. Default: 0.0001\n
    --no-regress            Do not perform the regression analysis and simply\n
                            output all PRS.\n
    --score                 Method to handle missing genotypes. By default, \n
                            final scores are averages of valid per-allele \n
                            scores with missing genotypes contribute an amount\n
                            proportional to imputed allele frequency. To throw\n
                            out missing observations instead (decreasing the\n
                            denominator in the final average when this happens),\n
                            use the 'no_mean_imputation' modifier. Alternatively,\n
                            you can use the 'center' modifier to shift all scores\n
                            to mean zero. \n
    --upper         | -u    The final p-value threshold. Default: 0.5\n
\nPRSlice:\n
    --prslice               Perform PRSlice where the whole genome is first cut\n
                            into bin size specified by this option. PRSice will\n
                            then be performed on each bin. Bins are then sorted\n
                            according to the their R2. PRSice is then performed\n
                            again to find the best bin combination.\n
                            This cannot be performed together with PRSet
\nTarget File:\n
    --binary-target         Indicate whether the target phenotype\n
                            is binary or not. Either T or F should be\n
                            provided where T represent a binary phenotype.\n
                            For multiple phenotypes, the input should be\n
                            separated by comma without space. Default: T\n
    --keep                  File containing the sample(s) to be extracted from\n
                            the target file. First column should be FID and\n
                            the second column should be IID. If --ignore-fid is\n
                            set, first column should be IID\n
                            Mutually exclusive from --remove\n
    --pheno-file    | -f    Phenotype file containing the phenotype(s).\n
                            First column must be FID of the samples and\n
                            the second column must be IID of the samples.\n
                            When --ignore-fid is set, first column must\n
                            be the IID of the samples.\n
                            Must contain a header if --pheno-col is\n
                            specified\n
    --pheno-col             Headers of phenotypes to be included from the\n
                            phenotype file\n
    --prevalence    | -k    Prevalence of all binary trait. If provided\n
    --remove                will adjust the ascertainment bias of the R2.\n
                            Note that when multiple binary trait is found,\n
                            you must provide prevalence information for\n
                            all of them.\n
    --target        | -t    Target genotype file. Currently support\n
                            both BGEN and binary PLINK format. For \n
                            multiple chromosome input, simply substitute\n
                            the chromosome number with #. PRSice will\n
                            automatically replace # with 1-22\n
    --type                  File type of the target file. Support bed \n
                            (binary plink) and bgen format. Default: bed\n
\nMisc:\n
    --all                   Output PRS for ALL threshold. WARNING: This\n
                            will generate a huge file\n
    --ignore-fid            Ignore FID for all input. When this is set,\n
                            first column of most file will be assume to\n
                            be IID instead of FID\n
    --out           | -o    Prefix for all file output\n
    --perm                  Number of permutation to perform. This will\n
                            generate the empirical p-value for the BEST\n
                            threshold\n
    --seed          | -s    Seed used for permutation. If not provided,\n
    --print-snp             system time will be used as seed. When same\n
                            seed and same input is provided, same result\n
                            should be generated\n
    --thread        | -n    Number of thread use\n
    --help          | -h    Display this help message\n";

# Library handling --------------------------------------------------------

if(!exists('startsWith', mode='function')){
  startsWith <- function(x, prefix){
    return(substring(x, 1, nchar(prefix)) == prefix)
  }
}

libraries <- c("ggplot2", "data.table", "optparse", "methods", "tools")
found = FALSE
argv = commandArgs(trailingOnly = TRUE)
dir_loc = grep("--dir", argv)
if (length(dir_loc) != 0) {
  dir_loc = dir_loc + 1
  found = TRUE
}

# INSTALL_PACKAGE: Functions for automatically install all required packages
InstalledPackage <- function(package)
{
  available <- suppressMessages(
                  suppressWarnings( 
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
        writeLines(paste("Trying to install ", package, " in ", dir, "/lib", sep = ""))
      }
      suppressMessages(suppressWarnings( 
        install.packages( package, lib = paste(dir, "/lib", sep = ""), repos = "http://cran.rstudio.com/")
      ))
    }
    if (!InstalledPackage(package))
      return(FALSE)
  }
  return(TRUE)
}

for (library in libraries)
{
  if (found)
  {
    if (!UsePackage(library, argv[dir_loc]))
    {
      stop("Error: ", library, " cannot be load nor install!")
    }
  } else{
    if (!UsePackage(library, "."))
    {
      stop("Error: ", library, " cannot be load nor install!")
    }
    
  }
}


# Command line arguments --------------------------------------------------

# We don't type the help message here. Will just directly use the usage information from c++
# See the Help Messages section for more information
option_list <- list(
  make_option(c("-b", "--base"), type = "character"),
  make_option(c("-B", "--bed"), type = "character"),
  make_option(c("-c", "--cov-header"), type = "character", dest="cov_header"),
  make_option(c("-C", "--cov-file"), type = "character", dest="cov_file"),
  make_option(c("-f", "--pheno-file"), type = "character", dest="pheno_file"),
  make_option(c("-g", "--gtf"), type = "character"),
  make_option(c("-i", "--interval"), type = "numeric"),
  make_option(c("-k", "--prevalence"), type = "numeric"),
  make_option(c("-l", "--lower"), type = "numeric"),
  make_option(c("-L", "--ld"), type = "character"),
  make_option(c("-m", "--msigdb"), type = "character"),
  make_option(c("-n", "--thread"), type = "numeric"),
  make_option(c("-o", "--out"), type = "character", default="PRSice"),
  make_option(c("-p", "--pvalue"), type = "character"),
  make_option(c("-s", "--seed"), type = "numeric"),
  make_option(c("-t", "--target"), type = "character"),
  make_option(c("-u", "--upper"), type = "numeric"),
  make_option(c("--all"), action = "store_true"),
  make_option(c("--beta"), action = "store_true"),
  make_option(c("--fastscore"), action = "store_true"),
  make_option(c("--full"), action = "store_true"),
  make_option(c("--ignore-fid"), action = "store_true", dest="ignore_fid"),
  make_option(c("--index"), action = "store_true"),
  make_option(c("--logit-perm"), action="store_true", dest="logit_perm"),
  make_option(c("--no-clump"), action = "store_true", dest="no_clump"),
  make_option(c("--no-regress"), action = "store_true", dest="no_regress"),
  make_option(c("--no-x"), action = "store_true", dest="no_x"),
  make_option(c("--no-y"), action = "store_true", dest="no_y"),
  make_option(c("--no-xy"), action = "store_true", dest="no_xy"),
  make_option(c("--no-mt"), action = "store_true", dest="no_mt"),
  make_option(c("--print-snp"), action = "store_true", dest="print_snp"),
  make_option(c("--A1"), type = "character"),
  make_option(c("--A2"), type = "character"),
  make_option(c("--bar-levels"), type = "character", dest="bar_levels", default="0.001,0.05,0.1,0.2,0.3,0.4,0.5"),
  make_option(c("--binary-target"), type = "character",dest="binary_target"),
  make_option(c("--bp"), type = "character"),
  make_option(c("--chr"), type = "character"),
  make_option(c("--clump-kb"), type = "character",dest="clump_kb"),
  make_option(c("--clump-p"), type = "numeric",dest="clump_p"),
  make_option(c("--clump-r2"), type = "numeric",dest="clump_r2"),
  make_option(c("--feature"), type = "character"),
  make_option(c("--keep"), type = "character"),
  make_option(c("--ld-type"), type = "character",dest="ld_type"),
  make_option(c("--num-auto"), type = "numeric", dest="num_auto"),
  make_option(c("--perm"), type = "numeric"),
  make_option(c("--pheno-col"), type = "character", dest="pheno_col"),
  make_option(c("--proxy"), type = "numeric"),
  make_option(c("--prslice"), type = "numeric"),
  make_option(c("--remove"), type = "character"),
  make_option(c("--score"), type = "character"),
  make_option(c("--se"), type = "character"),
  make_option(c("--snp"), type = "character"),
  make_option(c("--stat"), type = "character"),
  make_option(c("--type"), type = "character"),
  #R Specified options
  make_option(c("--plot"), action = "store_true"),
  make_option(c("--quantile", "-q"), type = "numeric"),
  make_option(c("--quant-extract", "-e"), type = "character", dest="quant_extract"),
  make_option("--quant-ref", type = "numeric", dest="quant_ref"),
  make_option("--scatter-r2", action = "store_true", default = F, dest="scatter_r2" ),
  make_option( "--bar-col-p", action = "store_true", default = F, dest="bar_col_p" ),
  make_option( "--bar-col-low", type = "character", default = "dodgerblue", dest="bar_col_low" ),
  make_option( "--bar-col-high", type = "character", default = "firebrick", dest="bar_col_high" ),
  make_option( "--bar-palatte", type="character", default = "YlOrRd", dest="bar_palatte" ),
  make_option("--prsice", type = "character"),
  make_option("--dir", type = "character")
  )


capture <- commandArgs(trailingOnly=TRUE)
help <- (sum(c("--help", "-h") %in%capture)>=1)
if(help){
  cat(help_message);
  quit();
}
argv <- parse_args(OptionParser(option_list = option_list))

not_cpp <-
  c(
    "help",
    "plot",
    "quantile",
    "quant-extract",
    "intermediate",
    "quant-ref",
    "scatter-r2",
    "bar-col-p",
    "bar-col-low",
    "bar-col-high",
    "bar-palatte",
    "prsice",
    "dir"
  )


# Check help messages --------------------------------------------------

provided <- function(name, argv) {
  return(name %in% names(argv))
}



# CALL_PRSICE: Call the cpp PRSice if required
# To ensure the excutable is set correctly
# This is also one of the reason why window doesn't work. I don't know if we can handle the \
# WINDOW PEOPLE
if (provided("prsice", argv)) {
  if (!startsWith(argv$prsice, "/") &&
      !startsWith(argv$prsice, ".")) {
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
command = ""
argv_c = argv
names(argv_c) = gsub("_", "-", names(argv))
flags <- c("all", "beta", "ignore-fid", "index", "no-clump", "no-regress", "no-x", "no-y", "no-xy", "no-mt",
           "fastscore", "print-snp")
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
  } else{
    stop("Cannot run PRSice without the PRSice binary file")
  }
}


# Plottings ---------------------------------------------------------------


multiplot <-
  function(...,
           plotlist = NULL,
           file,
           cols = 1,
           layout = NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                       ncol = cols,
                       nrow = ceiling(numPlots / cols))
    }
    
    if (numPlots == 1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]],
              vp = viewport(
                layout.pos.row = matchidx$row,
                layout.pos.col = matchidx$col
              ))
      }
    }
  }

# PLOTTING: Here contains all the function for plotting
# quantile_plot: plotting the quantile plots
quantile_plot <-
  function(PRS, PRS.best, pheno, prefix, argv, binary) {
    writeLines("Plotting the quantile plot")
    num_cov <- ncol(pheno)-2
    if(!provided("ignore_fid", argv)){
      num_cov <- num_cov-1
    }
    extract = NULL
    if (provided("quant_extract", argv)) {
      extract = fread(argv$quant_extract,
                      header = F,
                      data.table = F)
    }
    num_quant <- argv$quantile
    # Need to check if we have less pehnotypes than quantile
    colnames(PRS.best)[ncol(PRS.best)] <- "PRS"
    pheno.include <- NULL #Because we always name the phenotype as pheno, it will never be PRS
    
    if(provided("ignore_fid", argv)){
      pheno.merge<-merge(PRS.best, pheno, by="IID")
    }else{
      pheno.merge <- merge(PRS.best, pheno, by=c("FID", "IID"))
    }
    if (length(unique(pheno.merge$PRS)) < num_quant) {
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
    quants <-
      as.numeric(cut(
        pheno.merge$PRS,
        breaks = unique(quantile(pheno.merge$PRS, probs = seq(0, 1, 1 / num_quant))),
        include.lowest = T
      ))
    
    if (anyDuplicated(quantile(pheno.merge[, 4], probs = seq(0, 1, 1 / num_quant)))) {
      writeLines(paste(
        "Duplicate quantiles formed. Will use less quantiles: ",
        length(unique(quants)),
        sep = ""
      ))
      
    }
    num_quant = sum(!is.na(unique(quants)))
    if (!is.null(extract)) {
      extract_ID <- paste(extract$V1,extract$V2,sep="_")
      best_ID <- paste(pheno.merge$FID,pheno.merge$IID,sep="_")
      quants[best_ID %in% extract_ID] <- num_quant + 1 # We only matched based on the IID here
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
    
    if (num_cov > 0) {
      pheno.merge <-
        pheno.merge[, c("Pheno", "quantile", paste("Cov", 1:num_cov))]
    } else{
      pheno.merge <- pheno.merge[, c("Pheno", "quantile")]
    }
    
    
    family <- gaussian
    if (binary) {
      family <- binomial
    }
    reg <- summary(glm(Pheno ~ ., family, data = pheno.merge))
    coef.quantiles <- reg$coefficients[1:num_quant, 1]
    ci.quantiles.u <-
      reg$coefficients[1:num_quant, 1] + (1.96 * reg$coefficients[1:num_quant, 2])
    ci.quantiles.l <-
      reg$coefficients[1:num_quant, 1] - (1.96 * reg$coefficients[1:num_quant, 2])
    coef.quantiles[1] <- 0
    ci.quantiles.u[1] <- 0
    ci.quantiles.l[1] <- 0
    quantiles.for.table <-
      c(quant.ref, seq(1, num_quant, 1)[-quant.ref])
    quantiles.df <-
      data.frame(coef.quantiles,
                 ci.quantiles.u,
                 ci.quantiles.l,
                 quantiles.for.table)
    names(quantiles.df) <- c("Coef", "CI.U", "CI.L", "DEC")
    quantiles.df$Group = 0
    
    if (!is.null(extract)) {
      quantiles.df$Group[max(quantiles.df$DEC)] = 1
    }
    quantiles.df$Group <-
      factor(quantiles.df$Group, levels = c(0, 1))
    quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
    quantiles.plot <- ggplot(quantiles.df) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5)
      ) +
      ylab("Change in Phenotype given score in quantiles") +
      xlab("Quantiles for Polygenic Score") +
      scale_x_continuous(breaks = seq(0, num_quant, 1)) +
      theme(axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"))
    if (is.null(extract)) {
      quantiles.plot <-
        quantiles.plot + geom_point(aes(x = DEC, y = Coef),
                                    colour = "royalblue2",
                                    size = 4) +
        geom_pointrange(aes(
          ymin = CI.L,
          ymax = CI.U,
          y = Coef,
          x = DEC
        ),
        colour = "royalblue2",
        size = 0.9)
    } else{
      quantiles.plot <-
        quantiles.plot + geom_point(aes(x = DEC,
                                        y = Coef,
                                        color = Group), size = 4) +
        geom_pointrange(aes(
          ymin = CI.L,
          ymax = CI.U,
          y = Coef,
          x = DEC,
          color = Group
        ),
        size = 0.9) +
        scale_colour_manual(values = c("#0072B2", "#D55E00"))
    }
    ggsave(
      paste(prefix, "QUANTILES_PLOT.png", sep = "_"),
      width = 7,
      height = 7
    )
  }

high_res_plot <- function(PRS, prefix, argv) {
  # we will always include the best threshold
  writeLines("Plotting the high resolution plot")
  
  barchart.levels <-
    c(strsplit(argv$bar_level, split = ",")[[1]], PRS$Threshold[which.max(PRS$R2)])
  barchart.levels <-
    as.numeric(as.character(sort(unique(barchart.levels), decreasing = F)))
  # As the C++ program will skip thresholds, we need to artificially add the correct threshold information
  PRS.ori = PRS
  threshold <- as.numeric(as.character(PRS.ori$Threshold))
  for (i in 1:length(barchart.levels)) {
    if (sum(barchart.levels[i] - threshold > 0) > 0) {
      target <- max(threshold[barchart.levels[i] - threshold >= 0])
      temp <- PRS.ori[threshold == target, ]
      temp$Threshold <- barchart.levels[i]
      PRS = rbind(PRS, temp)
      
    } else{
      target <-
        (threshold[which(abs(threshold - barchart.levels[i]) == min(abs(threshold - barchart.levels[i])))])
      temp <- PRS.ori[threshold == target, ]
      temp$Threshold <- barchart.levels[i]
      PRS = rbind(PRS, temp)
      
    }
  }
  PRS = unique(PRS)
  # Need to also plot the barchart level stuff with green
  ggfig.points <- NULL
  if (argv$scatter_r2) {
    ggfig.points <- ggplot(data = PRS, aes(x = Threshold, y = R2)) +
      geom_line(aes(Threshold,  R2),
                colour = "green",
                data = PRS[with(PRS, Threshold %in% barchart.levels) , ]) +
      geom_hline(yintercept = max(PRS$R2), colour = "red") +
      ylab(expression(paste("PRS model fit:  ", R ^ 2, sep = " ")))
  } else{
    ggfig.points <-
      ggplot(data = PRS, aes(x = Threshold, y = -log10(P))) +
      geom_line(aes(Threshold, -log10(P)),
                colour = "green",
                data = PRS[with(PRS, Threshold %in% barchart.levels) , ]) +
      geom_hline(yintercept = max(-log10(PRS$P)), colour = "red") +
      ylab(bquote(PRS ~ model ~ fit:~ italic(P) - value ~ (-log[10])))
  }
  ggfig.points <- ggfig.points + geom_point() + geom_line() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", size = 0.5),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black")
    ) +
    xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T])))
  
  ggsave(
    paste(prefix, "_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""),
    width = 7,
    height = 7
  )
}

bar_plot <- function(PRS, prefix, argv) {
  writeLines("Plotting Bar Plot")
  barchart.levels <-
    c(strsplit(argv$bar_level, split = ",")[[1]], PRS$Threshold[which.max(PRS$R2)])
  barchart.levels <-
    as.numeric(as.character(sort(unique(barchart.levels), decreasing = F)))
  threshold <- as.numeric(as.character(PRS$Threshold))
  
  PRS.ori = PRS
  threshold <- as.numeric(as.character(PRS.ori$Threshold))
  for (i in 1:length(barchart.levels)) {
    if (sum(barchart.levels[i] - threshold > 0) > 0) {
      target <- max(threshold[barchart.levels[i] - threshold >= 0])
      temp <- PRS.ori[threshold == target, ]
      temp$Threshold <- barchart.levels[i]
      PRS = rbind(PRS, temp)
      
    } else{
      target <-
        (threshold[which(abs(threshold - barchart.levels[i]) == min(abs(threshold - barchart.levels[i])))])
      temp <- PRS.ori[threshold == target, ]
      temp$Threshold <- barchart.levels[i]
      PRS = rbind(PRS, temp)
      
    }
  }
  PRS <- unique(PRS[order(PRS$Threshold), ])
  # As the C++ program will skip thresholds, we need to artificially add the correct threshold information
  output <- PRS[PRS$Threshold %in% barchart.levels,]
  output$print.p[round(output$P, digits = 3) != 0] <-
    round(output$P[round(output$P, digits = 3) != 0], digits = 3)
  output$print.p[round(output$P, digits = 3) == 0] <-
    format(output$P[round(output$P, digits = 3) == 0], digits = 2)
  output$sign <- sign(output$Coefficient)
  output$print.p <- sub("e", "*x*10^", output$print.p)
  ggfig.plot <- ggplot(data = output)
  if (argv$bar_col_p) {
    ggfig.plot <-
      ggfig.plot + geom_bar(aes(
        x = factor(Threshold),
        y = R2,
        fill = factor(Threshold)
      ), stat = "identity") +
      scale_fill_brewer(palette = argv$palatte,
                        name = expression(italic(P) - value ~ threshold))
  }
  if (!argv$bar_col_p) {
    ggfig.plot <-
      ggfig.plot + geom_bar(aes(
        x = factor(Threshold),
        y = R2,
        fill = -log10(P)
      ), stat = "identity") +
      scale_fill_gradient(
        low = argv$bar_col_low,
        high = argv$bar_col_high,
        name = bquote(atop(-log[10] ~ model, italic(P) - value),)
      )
  }
  
  ggfig.plot <-
    ggfig.plot + geom_text(
      aes(
        x = factor(Threshold),
        y = R2,
        label = paste(print.p)
      ),
      vjust = -1.5,
      hjust = 0,
      angle = 45,
      cex = 2.8,
      parse = T
    ) +
    scale_y_continuous(limits = c(0, max(output$R2) * 1.25)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", size = 0.5) ,
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black")
    ) +
    xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
    ylab(expression(paste("PRS model fit:  ", R ^ 2)))
  ggsave(
    paste(prefix, "_BARPLOT_", Sys.Date(), ".png", sep = ""),
    width = 7,
    height = 7
  )
}

# run_plot: The function used for calling different plotting functions
run_plot <- function(prefix, argv, pheno_matrix, binary) {
  writeLines("")
  #writeLines(prefix)
  PRS <- fread(paste(prefix, ".prsice", sep = ""), header = T, data.table = F)
  PRS.best <- fread(paste(prefix, ".best", sep = ""), header = T, data.table = F)
  
  # start from here, we need to organize all the file accordingly so that the individual actually match up with each other
  # Good thing is, only quantile plot really needs the cov and phenotype information
  if (provided("quantile", argv) && argv$quantile > 0) {
    # Main purpose, match up with PRS.best as that is the input
    PRS.best.reduce <- subset(PRS.best, PRS.best$Included=="Y")
    # Need to plot the quantile plot (Remember to remove the iid when performing the regression)
    quantile_plot(PRS, PRS.best.reduce, pheno_matrix, prefix, argv, binary)
  }
  # Now perform the barplotting
  if (!provided("fastscore", argv) || !argv$fastscore) {
    high_res_plot(PRS, prefix, argv)
  }
  bar_plot(PRS, prefix, argv)
}




# Process file names for plotting------------------------------------------------------
if(provided("no_regress", argv)){
    quit("yes");
}
ignore_fid <- provided("ignore_fid", argv)
extract_matrix <- function(x,y){
  z=which(x==y); 
  if(length(z)==0){
    return(NA)
  }else{
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

if(!provided("pheno_file", argv) && !provided("binary_target", argv)){
  argv$binary_target="T";
}else{
  if(!provided("pheno_col",argv) && provided("binary_target", argv) 
     && length(argv$binary_target)==1){
    #This is ok
    test <- "ok"
  }else if(!provided("pheno_col", argv) && !provided("binary_target",argv)){
    argv$binary_target="T";
  }else if(provided("pheno_col", argv) && !provided("binary_target", argv)
           && length(argv$pheno_col)<=1){
    argv$binary_target = "T"
  }
  else if(provided("pheno_col", argv) && provided("binary_target", argv)
          && length(argv$pheno_col)!=length(argv$binary_target)){
    stop("ERROR: Number of target phenotypes doesn't match information of binary target! 
    You must indicate whether the phenotype is binary using --binary-target")
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
  }else{
    header <- read.table(argv$pheno_file, nrows = 1, header = TRUE)
    # This will automatically filter out un-used phenos
    if (length(binary_target) != length(phenos)) {
      message <-"Number of binray target should match number of phenotype provided!"
      message = paste( message,
            "There are ", length(binary_target), " binary target information and ",
            length(phenos), "phenotypes", sep = "")
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


# Now we have the correct header of phenos, bases and binary_target information
# Need to get the covariates
covariance = NULL
if (provided("cov_file", argv)) {
  covariance = fread(argv$cov_file, data.table = F, header = T)
  if (provided("cov_header", argv)) {
    c = strsplit(argv$cov_header, split = ",")[[1]]
    selected = colnames(covariance)%in%c
    if(!ignore_fid){
      selected[2] <- TRUE #When ignore_fid isn't provided, then we need to also include the FID information
    }
    selected[1] <- TRUE # we always want the IID information, which should be the first column
    covariance <- covariance[, selected]
  }
  if(ignore_fid){
    colnames(covariance)<-c("IID", paste("Cov",1:(ncol(covariance)-1)))
  }else{
    colnames(covariance)<-c("FID","IID", paste("Cov",1:(ncol(covariance)-2)))
  }
}

# we no longer have those complication
prefix <- argv$out

num_region = nrow(read.table(paste(prefix, "region", sep="."), header=T));

region = NULL
# we need the fam file under all situation
fam = fread(
paste(argv$target, ".fam", sep = ""), data.table = F, header = F )
colnames(fam)[1:2] = c("FID", "IID")
fam_id <- paste(fam$FID, fam$IID, sep="_")
match_cov = NULL
if (!is.null(covariance)) {
  if(ignore_fid){
    match_cov <- covariance[covariance$IID %in% fam$IID,]
    match_cov <- match_cov[match(fam$IID, match_cov$IID),] # match the ordering
    match_cov <- match_cov[ !is.na(apply(match_cov[,-1],1,sum)),]
  }else{
    cov_ID <- paste(covariance$FID, covariance$IID, sep="_")
    match_cov = covariance[cov_ID %in% fam_id,]
    cov_ID <- paste(match_cov$FID, match_cov$IID, sep="_") #updated cov_ID
    match_cov = match_cov[match(fam_id, cov_ID),] # match the ordering
    match_cov <- match_cov[ !is.na(apply(match_cov[,-c(1:2)],1,sum)),] #remove all NA
  }
}
#Now match_cov contain all samples with valid covariates
if (!is.null(phenos)) {
  pheno_file = fread(argv$pheno_file, header = T, data.table = F)
  if(ignore_fid){
    colnames(pheno_file)[1] <- "IID"
  }else{
    colnames(pheno_file)[1:2] <- c("FID", "IID")
  }
  if(ignore_fid){
    #Only care about the IID
    pheno_file <- pheno_file[pheno_file$IID %in% fam$IID, ] # This will select only samples found within the fam file
    pheno_file <- pheno_file[match(match_id, pheno_file$IID),] # match the ordering
  }else{
    pheno_id <- paste(pheno_file$FID, pheno_file$IID, sep="_")
    pheno_file <- pheno_file[pheno_id %in% fam_id,]
    pheno_id <- paste(pheno_file$FID, pheno_file$IID, sep="_") #Updated pheno_id
    pheno_file <- pheno_file[match(fam_id, pheno_id),]
  }
 
  id = 1
  phenos.index <- unlist(apply(as.matrix(phenos), 1, extract_matrix ,colnames(pheno_file)))
        
  for (p in 1:length(phenos.index)) {
    if(!is.na(phenos.index[p])){
      cur_prefix = paste(prefix, phenos[p], region, sep = ".")
      if(num_region==1){
        cur_prefix = paste(prefix, phenos[p], sep = ".")
      }
      # Give run_plot a ready to use matrix
      cur_pheno = NULL
      if(ignore_fid){
        cur_pheno = data.frame(IID=pheno_file[,1], Pheno=pheno_file[phenos.index[p]])
      }else{
        cur_pheno = data.frame(FID=pheno_file[,1], IID=pheno_file[,2], Pheno=pheno_file[phenos.index[p]])
      }
      if(binary_target[id]){
        # Update the cur_pheno
        cur_pheno$Pheno = suppressWarnings(as.numeric(as.character(cur_pheno$Pheno)))
        cur_pheno$Pheno[cur_pheno$Pheno > 2] = NA
        cur_pheno$Pheno[cur_pheno$Pheno < 0] = NA
        if(max(cur_pheno$Pheno, na.rm=T)==2 && min(cur_pheno$Pheno, na.rm=T)==0){
          stop("Invalid case control formating. Either use 0/1 or 1/2 coding")
        }else if(max(cur_pheno$Pheno, na.rm=T)==2){
          cur_pheno$Pheno = cur_pheno$Pheno-1
        }
      }
      cur_pheno.clean = cur_pheno[!is.na(cur_pheno$Pheno),]
      fam.final = cur_pheno.clean
      if(!is.null(match_cov)){
        if(ignore_fid){
          fam.final = merge(cur_pheno.clean, match_cov, by="IID");
        }else{
          fam.final = merge(cur_pheno.clean, match_cov, by=c("FID", "IID"));
        }
      }
      fam.ok = fam[fam$IID%in%fam.final$IID,]
      fam.final = fam.final[match(fam.ok$IID, fam.final$IID),]
      run_plot(cur_prefix,argv,fam.final,binary_target[id])
    }else{
      writeLines(paste(phenos[p],"not found in the phenotype file. It will be ignored"))
    }
    id = id + 1
  }
} else{
  # No phenotype headers
  cur_prefix = paste(prefix, region, sep = ".")
  if(num_region==1){
    cur_prefix = paste(prefix, sep = ".")
  }
  pheno=NULL
  if (provided("pheno_file", argv)) {
    pheno = fread(paste(argv$pheno_file), data.table = F, header = F) 
    if(ignore_fid){
      colnames(pheno)[1] <- "IID"
    }else{
      colnames(pheno)[1:3] <- c("FID", "IID", "V2") #Otherwise this is V3
    }
    #Unless someone being stupid and name their sample's FID and IID as the header, this should be fine
  }
  # Give run_plot a ready to use matrix
  fam.clean = data.frame(FID=fam$FID, IID=fam$IID, Pheno=fam$V6);
  if(!is.null(pheno))
  {
    if(ignore_fid){
      fam.clean = data.frame(IID=pheno$IID, Pheno=pheno$V2)
    }else{
      fam.clean = data.frame(FID=pheno$FID, IID=pheno$IID, Pheno=pheno$V2)
    }
  }
  if(binary_target[1]){
    # Update the cur_pheno
    fam.clean$Pheno = suppressWarnings(as.numeric(as.character(fam.clean$Pheno)))
    fam.clean$Pheno[fam.clean$Pheno > 2] = NA
    fam.clean$Pheno[fam.clean$Pheno < 0] = NA
    if(max(fam.clean$Pheno, na.rm=T)==2 && min(fam.clean$Pheno, na.rm=T)==0){
      stop("Invalid case control formating. Either use 0/1 or 1/2 coding")
    }else if(max(fam.clean$Pheno, na.rm=T)==2){
      fam.clean$Pheno = fam.clean$Pheno-1
    }
  }
  cur_pheno.clean = fam.clean[!is.na(fam.clean$Pheno),]
  fam.final = cur_pheno.clean
  if(!is.null(match_cov)){
    if(ignore_fid){
      fam.final <- merge(cur_pheno.clean, match_cov, by="IID")
    }else{
      fam.final <- merge(cur_pheno.clean, match_cov, by=c("FID", "IID"))
    }
  }
  fam.ok = fam[fam$IID%in%fam.final$IID,]
  fam.final = fam.final[match(fam.ok$IID, fam.final$IID),]
  run_plot(cur_prefix, argv, fam.final, binary_target[1])
}
