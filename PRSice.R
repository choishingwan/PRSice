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

# Environment stuff: This will allow us to locate the cpp file correctly


# Library handling --------------------------------------------------------

if(!exists('startsWith', mode='function')){
  startsWith <- function(x, prefix){
    return(substring(x, 1, nchar(prefix)) == prefix)
  }
}

libraries <-
  c("ggplot2", "data.table", "optparse", "methods", "tools")
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
  available <-
    suppressMessages(suppressWarnings(
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
        writeLines(paste("Trying to install ",
                         package,
                         " in ",
                         dir,
                         "/lib",
                         sep = ""))
      }
      suppressMessages(suppressWarnings(
        install.packages(
          package,
          lib = paste(dir, "/lib", sep = ""),
          repos = "http://cran.rstudio.com/"
        )
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


option_list <- list(
  make_option(c("-b", "--base"), type = "character", help = "Base association files. User can provide multiple base files"),
  make_option(
    c("-t", "--target"),
    type = "character",
    help = "Plink binary file prefix for target files. Currently only
    support plink binary inputs. For multiple target phenotypes, user should use the --pheno_file option
    together with the pheno_col option. For multiple chromosome input, one  should substitute the chromosome
    number with #. For example, if you files are presented as genotype_chr1_test, genotype_chr2_test, genotype_chr3_test,
    then you can use: genotype_chr#_test. Please note that the substitute is based on your base file. So if
    your base file code chromosome with chr, e.g. chr1 chr2 etc, then in our example case, you should code your plink
    file as genotype_#_test"
  ),
  make_option(
    "--binary_target",
    type = "character",
    help = "Indicate whether the target sample has binary phenotype or not.
    For each phenotype, user need to provide either T or F where T means the phenotype is binary"
  ),
  
  make_option("--beta", help = "Indicate whether the test statistic is beta instead of OR. Must be of the same length as base"),
  make_option(
    c("--pheno_file", "-f"),
    type = "character",
    help = "Phenotype file containing the target phenotype(s).
    If provided, the fam file of the target is ignored. First column must be IID of the samples. Must contain a header
    if --pheno_col is specified"
  ),
  make_option("--pheno_col", type = "character", help = "Headers of pheenotypes from phenotype file"),
  make_option(
    c("--ld", "-L"),
    type = "character",
    help = "Plink binary file prefix for the reference file used for LD calculation.
    If not provided, will use the target genotype for the LD calculation. Can also use multiple chromosome
    plink file. Please see --target for more information."
  ),
  make_option(
    c("--covar_header", "-c"),
    type = "character",
    help = "Header of covariates. If not provided, will use all variable in the
    covariate file as the covarite."
  ),
  make_option(c("--covar_file", "-C"), type = "character", help = "Covarite file. Formate should be: ID Cov1 Cov2 ..."),
  make_option("--full", action = "store_true", help = "Also include the full model in the PRSice output"),
  make_option("--all", action = "store_true", help = "Output PRS for ALL threshold. Can only be used together with fastscore to avoid huge output files."),
  make_option(
    "--no_regress",
    action = "store_true",
    help = "Do not perform the regression analysis and simply output all PRS. Can only be used together
    with fastscore to avoid huge output files. If you must, you can modify bar_levels to obtain the fine scale PRS outputs"
  ),
  make_option(
    c("--out", "-o"),
    type = "character",
    help = "Prefix of all output.",
    default = "PRSice"
  ),
  make_option(
    c("--lower", "-l"),
    type = "numeric",
    help = "The starting p-value threshold. ",
    default = 0.0001
  ),
  make_option(
    c("--upper", "-u"),
    type = "numeric",
    help = "The final p-value threshold.",
    default = 0.5
  ),
  make_option(
    c("--interval", "-i"),
    type = "numeric",
    help = "The step size of the threshold.",
    default = 0.00005
  ),
  make_option("--fastscore", action = "store_true", help = "Calculate the minimum amount of threshold as required by the bar_level option"),
  make_option("--chr", type = "character", help = "Column header of Chromosome <Required>"),
  make_option("--A1", type = "character", help = "Column header of Reference Allele <Required>"),
  make_option("--A2", type = "character", help = "Column header of Alternaative Allele"),
  make_option("--stat", type = "character", help = "Column header of test statistic <Required>"),
  make_option("--snp", type = "character", help = "Column header of SNP id"),
  make_option("--bp", type = "character", help = "Column header of SNP location"),
  make_option("--se", type = "character", help = "Column header of Standard Error"),
  make_option(c("--pvalue", "-p"), help = "Column header of p-value <Required> "),
  make_option(
    "--index",
    action = "store_true",
    help = "Indicate all the above options are providing the INDEX of the corresponding column.
    (Index should be 0-based). Useful when your base file each have a different header but the column
    index remains the same"
  ),
  make_option(
    "--clump_p",
    type = "numeric",
    help = "The p-value threshold use for clumping.",
    default = 1
  ),
  make_option(
    "--clump_r2",
    type = "numeric",
    help = "The R2 threshold for clumping. Please note that as we did not implement the
    maximum likelihood R2 calculation, the clumping result can differ slightly from plink.",
    default = 0.1
  ),
  make_option(
    "--clump_kb",
    type = "numeric",
    help = "The distance for clumping in kb.",
    default = 250
  ),
  make_option(c("--bed", "-B"), type = "character", help = "Bed file containing the selected regions. Name of bed file will be used as the region identifier."),
  make_option(c("--gtf", "-g"), type = "character", help = "GTF file containing gene boundaries. Required when --msigdb is set."),
  make_option(c("--msigdb", "-m"), type = "character", help = "MSIGDB file containing the pathway information require the gtf file."),
  make_option("--gen_bed", action = "store_true", help = "Generate bed file of gene regions from the gtf file."),
  make_option(
    "--proxy",
    type = "numeric",
    help = "Proxy threshold for index SNP to be considered as part of the region represented by the clumped SNPs.
    e.g. --proxy 0.8 means the index SNP will represent the region of any clumped SNPs that has a R2 >= 0.8 with it even if
    it is not physically within these regions"
  ),
  make_option("--feature", type="character", help="Features to be included from the gtf file. Default is exon, CDS, gene and protein_coding. If this parameter is provided, all default will be ignored."),
  make_option(
    "--prslice",
    type = "numeric",
    help = "Perform PRSlice where the whole genome is first cut into bin size specified by this option. PRSice
    will then be performed on each bin. Bins are then sorted according to the their R2. PRSice is then performed again
    to find the best bin combination. This cannot be performed together with PRSet"
  ),
  make_option(
    "--bar_levels",
    type = "character",
    help = "Level of barchart to be plotted. When fastscore is set, PRSice will
    only calculate the PRS for threshold within the bar level"
  ),
  make_option(
    c("--thread", "-T"),
    type = "numeric",
    help = "Number of thread use",
    default = 1
  ),
  make_option(
    "--perm", type="numeric", help="Number of permutation to perform. When this parameter is provided, permutation will be performed to obtain an empirical P-value. This will significantly increases the run time of PRSice."),
  make_option(
    "--c_help",
    action = "store_true",
    help = "Print the help message from the c++ program instead",
    default = F
  ),
  make_option(
    "--plot",
    action = "store_true",
    help = "Indicate only plotting is required",
    default = F
  ),
  make_option("--intermediate", type = "character", help = "Pefix of the intermediate files for plotting (e.g. ignore .prsice and .best). If not provided, will deduce the file prefix from --base and --target"),
  make_option(
    c("--quantile", "-q"),
    type = "numeric",
    help = "Number of quantiles to plot. 0 = Not producing the quantile plot",
    default = 0
  ),
  make_option(c("--quant_extract", "-e"), type = "character", help = "File containing sample ID to be plot on a separated quantile e.g. extra quantile containing only schizophrenia samples"),
  make_option(
    "--bar_level",
    type = "character",
    help = "barchar level used for plotting",
    default = "0.001,0.05,0.1,0.2,0.3,0.4,0.5"
  ),
  make_option("--quant_ref", type = "numeric", help = "Reference quantile for quantile plot"),
  make_option(
    "--scatter_r2",
    action = "store_true",
    help = "y-axis of the high resolution scatter plot should be R2",
    default = F
  ),
  make_option(
    "--bar_col_p",
    action = "store_true",
    help = "Change the colour of bar to p-value threshold instead of the association with phenotype",
    default = F
  ),
  make_option(
    "--bar_col_low",
    type = "character",
    help = "Colour of the poorest predicting threshold",
    default = "dodgerblue"
  ),
  make_option(
    "--bar_col_high",
    type = "character",
    help = "Colour of the most predicting threshold",
    default = "firebrick"
  ),
  make_option(
    "--bar_palatte",
    type="character",
    help ="Colour palatte to be used for bar plotting when --bar_col_p is set",
    default = "YlOrRd"
  ),
  make_option("--prsice", type = "character", help = "Location of the PRSice binary"),
  make_option("--dir", type = "character", help = "Location to install ggplot. Only require if ggplot is not installed")
  )

argv <- parse_args(OptionParser(option_list = option_list))
# stop()

# argv = commandArgs(trailingOnly = TRUE)
#help=(sum(c("--help", "-h") %in%argv)>=1)
#if(help){
#  print.arg.parser(p)
#  quit();
#}
#argv <- parse_args(p)
not_cpp <-
  c(
    "help",
    "c_help",
    "plot",
    "quantile",
    "quant_extract",
    "intermediate",
    "quant_ref",
    "scatter_r2",
    "bar_col_p",
    "bar_col_low",
    "bar_col_high",
    "bar_palatte",
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
if (provided("prsice", argv)) {
  if (!startsWith(argv$prsice, "/") &&
      !startsWith(argv$prsice, ".")) {
    argv$prsice = paste("./", argv$prsice, sep = "")
  }
}

if (argv$c_help) {
  if (!provided("prsice", argv)) {
    stop("Cannot use c_help without specifying the location of the PRSice binary!")
    
  }
  system(paste(argv$prsice, " --help", sep = ""))
  quit()
  
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
if (!argv$plot) {
  for (i in names(argv)) {
    # only need special processing for flags and specific inputs
    if (i == "index" ||
        i == "gen_bed" ||
        i == "fastscore" ||
        i == "full" || i == "all" || i == "no_regress") {
      if (argv[[i]])
        command = paste(command, " --", i, sep = "")
    } else if (i %in% not_cpp) {
      # ignore
    } else{
      temp = add_command(argv[[i]])
      if (!is.na(temp)) {
        command = paste(command, " --", i, " ", temp, sep = "")
      }
    }
  }
  if (nchar(command) == 0) {
    print.arg.parser(p)
    quit()
  }
  if (provided("prsice", argv)) {
    ret <- system2(argv$prsice,
                   command,
                   stdout = TRUE,
                   stderr = TRUE)
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
    extract = NULL
    if (provided("quant_extract", argv)) {
      extract = fread(argv$quant_extract,
                      header = F,
                      data.table = F)
    }
    num_quant <- argv$quantile
    # Need to check if we have less pehnotypes than quantile
    if (length(unique(PRS.best[, 2])) < num_quant) {
      writeLines(
        paste(
          "WARNING: There are only ",
          length(unique(PRS.best[, 2])),
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
        PRS.best[, 2],
        breaks = unique(quantile(PRS.best[, 2], probs = seq(0, 1, 1 / num_quant))),
        include.lowest = T
      ))
    
    if (anyDuplicated(quantile(PRS.best[, 2], probs = seq(0, 1, 1 / num_quant)))) {
      writeLines(paste(
        "Duplicate quantiles formed. Will use less quantiles: ",
        length(unique(quants)),
        sep = ""
      ))
      
    }
    num_quant = sum(!is.na(unique(quants)))
    if (!is.null(extract)) {
      quants[PRS.best[, 1] %in% extract$V2] <- num_quant + 1
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
    pheno$quantile <- quants
    
    if (ncol(pheno) > 2) {
      pheno <-
        pheno[, c(colnames(pheno)[1], "quantile", colnames(pheno)[2:(ncol(pheno) -
                                                                       1)])]
    } else{
      pheno <- pheno[, c(colnames(pheno)[1], "quantile")]
    }
    
    family <- gaussian
    if (binary) {
      family <- binomial
    }
    reg <- summary(glm(Pheno ~ ., family, data = pheno))
    
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
  writeLines(prefix)
  
  PRS <- fread(paste(prefix, ".prsice", sep = ""),
               header = T,
               data.table = F)
  PRS.best <-
    fread(paste(prefix, ".best", sep = ""),
          header = T,
          data.table = F)
  
  # start from here, we need to organize all the file accordingly so that the individual actually match up with each other
  # Good thing is, only quantile plot really needs the cov and phenotype information
  if (provided("quantile", argv) && argv$quantile > 0) {
    # Main purpose, match up with PRS.best as that is the input
    PRS.best = PRS.best[PRS.best$IID%in% pheno_matrix$IID, ]
    PRS.best = PRS.best[match(pheno_matrix$IID, PRS.best$IID),]
    # Need to plot the quantile plot (Remember to remove the iid)
    quantile_plot(PRS, PRS.best, pheno_matrix[,-1], prefix, argv, binary)
  }
  # Now perform the barplotting
  if (!provided("fastscore", argv) || !argv$fastscore) {
    high_res_plot(PRS, prefix, argv)
  }
  bar_plot(PRS, prefix, argv)
}




# Process file names for plotting------------------------------------------------------

# CALL PLOTTING FUNCTION: Process the input names and call the actual plotting function
if (provided("intermediate", argv)) {
  # File name provided
  # Need to know the beta and the phenotype of the target
  # So in theory, we still need argv$target
  # TODO: It is more complicated than this
  #run_plot(argv$intermediate, argv)
  stop("Currently not support intermediate. Use --plot instead")
} else{
  # we need to deduce the file names
  # Now we actually require one single string for the input, separated by ,
  # Get all the region information
  
  if (provided("base", argv) && !is.na(argv[["base"]])) {
    bases = strsplit(argv$base, split = ",")[[1]]
    # now check the target
    phenos = NULL
    if (!provided("target", argv)) {
      stop(
        "Target file name not found. You'll need to either provide the intermediate prefix of all the target name for plotting!"
      )
    }
    if (!provided("binary_target", argv)) {
      stop(
        "We do need to know if the target file is binary or not in order to decide whether if we will run logistic or linear regression."
      )
    }
    binary_target = strsplit(argv$binary_target, split = ",")[[1]]
    if (provided("pheno_col", argv)) {
      phenos = strsplit(argv$pheno_col, split = ",")[[1]]
      if (!provided("pheno_file", argv)) {
        phenos = NULL
        writeLines(
          strwrap(
            "WARNING: Cannot have multiple phenotypes if pheno_file is not provided. We will ignore the pheno_col option.",
            width = 80
          )
        )
      } else{
        header <- read.table(argv$pheno_file,
                             nrows = 1,
                             header = FALSE)
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
        binary_target = binary_target[phenos %in% header[1, ]]
        phenos = phenos[phenos %in% header[1, ]]
      }
    } else{
      if (length(binary_target) != 1) {
        stop("Too many binary target information. We only have one phenotype")
      }
    }
    # Now we have the correct header of phenos, bases and binary_target information
    # Need to get the covariates
    covariance = NULL
    if (provided("covar_file", argv)) {
      covariance = fread(argv$covar_file,
                         data.table = F,
                         header = T)
      if (provided("covar_header", argv)) {
        c = strsplit(argv$covar_header, split = ",")[[1]]
        selected = colnames(covariance)%in%c
        selected[1] = TRUE # we always want the IID information, which should be the first column
        covariance = covariance[, selected]
      }
    }
    for (b in bases) {
      base_name = file_path_sans_ext(basename(b))
      # region will always have at least one item -> Base
      prefix = paste(argv$out, base_name, sep = ".")
      region = "Base"
      
      if (!is.null(phenos)) {
        pheno_file = fread(argv$pheno_file,
                           header = T,
                           data.table = F)
        fam = fread(
          paste(argv$target, ".fam", sep = ""),
          data.table = F,
          header = F
        )
        pheno_file = pheno_file[pheno_file[, 1] %in% fam$V2, ] # This will select only samples found within the fam file
        pheno_file = pheno_file[match(fam$V2, pheno_file[,1]),] # match the ordering
        
        match_cov = NULL
        if (!is.null(covariance)) {
          match_cov = covariance[covariance[, 1] %in% fam$V2,]
          match_cov = match_cov[match(fam$V2, match_cov[,1]),] # match the ordering
        }
        id = 1
        phenos.index <- unlist(apply(as.matrix(phenos), 1, function(x,y){z=which(x==y); if(length(z)==0){return(NA)}else{return(z)}},colnames(pheno_file)))
        for (p in 1:length(phenos.index)) {
          if(!is.na(phenos.index[p])){
            cur_prefix = paste(prefix, phenos[p], region, sep = ".")
            # Give run_plot a ready to use matrix
            cur_pheno = data.frame(IID=pheno_file[,1], Pheno=pheno_file[phenos.index[p]])
            if(binary_target[id]){
              # Update the cur_pheno
              cur_pheno$Pheno = as.numeric(as.character(cur_pheno$Pheno))
              cur_pheno$Pheno[cur_pheno$Pheno > 2] = NA
              cur_pheno$Pheno[cur_pheno$Pheno < 0] = NA
              if(max(cur_pheno$Pheno, na.rm=T)==2 && min(cur_pheno$Pheno, na.rm=T)==0){
                stop("Invalid case control formating. Either use 0/1 or 1/2 coding")
              }else if(max(cur_pheno$Pheno, na.rm=T)==2){
                cur_pheno$Pheno = cur_pheno$Pheno-1
              }
            }
              
            cur_pheno.clean = cur_pheno[!is.na(cur_pheno$Pheno),]
            cov_matrix = match_cov[match_cov[,1]%in% cur_pheno.clean$IID & !is.na(apply(match_cov[,-1],1,sum)),]
            cur_pheno.clean = cur_pheno.clean[cur_pheno.clean$IID%in%cov_matrix[,1],]
            colnames(cov_matrix) = c("IID", paste("Cov",1:(ncol(cov_matrix)-1)))
            # Should be of same size now
            fam.final = merge(cur_pheno.clean, cov_matrix, by="IID");
            fam.ok = fam[fam$V2%in%fam.final$IID,]
            fam.final = fam.final[match(fam.ok$V2, fam.final$IID),]
            run_plot(cur_prefix,
                     argv,
                     fam.final,
                     binary_target[id])
          }else{
            writeLines(paste(phenos[p],"not found in the phenotype file. It will be ignored"))
          }
          id = id + 1
        }
      } else{
        cur_prefix = paste(prefix, region, sep = ".")
        
        fam = fread(
          paste(argv$target, ".fam", sep = ""),
          data.table = F,
          header = F
        )
        match_cov = NULL
        pheno=NULL
        if (!is.null(covariance)) {
          match_cov = covariance[covariance[, 1] %in% fam$V2,]
          match_cov = match_cov[match(fam$V2, match_cov[,1]), ]
        }
        if (provided("pheno_file", argv)) {
          pheno = fread(paste(argv$pheno_file),
                        data.table = F,
                        header = F)
        }
        # Doesn't matter if pheno file has a header or not

        # Give run_plot a ready to use matrix
        fam.clean = data.frame(IID=fam$V2, Pheno=fam$V6);
        if(!is.null(pheno))
        {
          fam.clean = data.frame(IID=pheno$V1, Pheno=pheno$V2)
        }
        if(binary_target[1]){
          # Update the cur_pheno
          fam.clean$Pheno = as.numeric(as.character(fam.clean$Pheno))
          fam.clean$Pheno[fam.clean$Pheno > 2] = NA
          fam.clean$Pheno[fam.clean$Pheno < 0] = NA
          if(max(fam.clean$Pheno, na.rm=T)==2 && min(fam.clean$Pheno, na.rm=T)==0){
            stop("Invalid case control formating. Either use 0/1 or 1/2 coding")
          }else if(max(fam.clean$Pheno, na.rm=T)==2){
            fam.clean$Pheno = fam.clean$Pheno-1
          }
        }
        cur_pheno.clean = fam.clean[!is.na(fam.clean$Pheno),]
        cov_matrix = match_cov[match_cov[,1]%in% cur_pheno.clean$IID & !is.na(apply(match_cov[,-1],1,sum)),]
        cur_pheno.clean = cur_pheno.clean[cur_pheno.clean$IID%in%cov_matrix[,1],]
        colnames(cov_matrix) = c("IID", paste("Cov",1:(ncol(cov_matrix)-1)))
        # Should be of same size now
        fam.final = merge(cur_pheno.clean, cov_matrix, by="IID");
        fam.ok = fam[fam$V2%in%fam.final$IID,]
        fam.final = fam.final[match(fam.ok$V2, fam.final$IID),]
        
        
        run_plot(cur_prefix, argv, fam.final, binary_target[1])
      }
    }
  } else{
    message <- "Base file name not found. "
    if (provided("plot", argv) && argv$plot) {
      message <-
        paste(
          message,
          "You'll need to either provide the intermediate prefix or all the target/base name for plotting"
        )
    }
    stop(message)
  }
}
