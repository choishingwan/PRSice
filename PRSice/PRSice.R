# Here is the guide to this protentially long R code
# To go to each section, just search for the corresponding header as stated here
# The code structure are as follow
#
# ARG
# The package argparser is included here to reduce one level of dependency
# This will allow us to parse the command line input elegantly, then decide
# where to install ggplot2. 
# END_ARG
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

libraries <- c( "ggplot2", "data.table", "optparse", "methods")
argv = commandArgs(trailingOnly = TRUE)
dir_loc = grep("--dir",argv)
if(length(dir_loc)!=0){
  dir_loc=dir_loc+1;
}
# INSTALL_PACKAGE: Functions for automatically install all required packages
InstalledPackage <- function(package)
{
  available <- suppressMessages(suppressWarnings(sapply(package, require, quietly = TRUE, character.only = TRUE, warn.conflicts = FALSE)))
  missing <- package[!available]
  if (length(missing) > 0) return(FALSE)
  return(TRUE)
}

CRANChoosen <- function()
{
  return(getOption("repos")["CRAN"] != "@CRAN@")
}

UsePackage <- function(package,dir)
{
  if(!InstalledPackage(package))
  {
    dir.create(file.path(dir,"lib"), showWarnings = FALSE)
    .libPaths(c(.libPaths(), paste(dir,"/lib",sep="")))
    if(!InstalledPackage(package)){
      if(is.na(dir)){
        print("WARNING: dir not provided, cannot install the required packages")
        return(FALSE);
      }else{
        print(paste("Trying to install ", package, " in ", dir, "/lib",sep=""))
      }
      suppressMessages(suppressWarnings(install.packages(package, lib=paste(dir, "/lib",sep=""),repos="http://cran.rstudio.com/")))
    }
    if(!InstalledPackage(package)) return(FALSE)
  }
  return(TRUE)
}

for(library in libraries)
{
  if(!UsePackage(library, argv[dir_loc]))
  {
    stop("Error: ", library, " cannot be load nor install!")
  }
}

option_list <- list(
  make_option(c("-b", "--base"), type="character", help="Base association files. User can provide multiple base files"),
  make_option(c("-t", "--target"), type="character", help="Plink binary file prefix for target files. Currently only 
                  support plink binary inputs. For multiple target phenotypes, user should use the --pheno_file option 
              together with the pheno_col option. For multiple chromosome input, one  should substitute the chromosome 
              number with #. For example, if you files are presented as genotype_chr1_test, genotype_chr2_test, genotype_chr3_test,
              then you can use: genotype_chr#_test. Please note that the substitute is based on your base file. So if 
              your base file code chromosome with chr, e.g. chr1 chr2 etc, then in our example case, you should code your plink 
              file as genotype_#_test"),
  make_option("--target_is_binary", type="character", help="Indicate whether the target sample has binary phenotype or not. 
                  For each phenotype, user need to provide either T or F where T means the phenotype is binary"),
  
  make_option("--beta", type="charactre", help="Indicate whether the test statistic is beta instead of OR. Must be of the same length as base"),
  make_option(c("--pheno_file","-f"), type="character", help="Phenotype file containing the target phenotype(s). 
                    If provided, the fam file of the target is ignored. First column must be IID of the samples. Must contain a header 
                    if --pheno_col is specified"),
  make_option("--pheno_col", type="character", help="Headers of pheenotypes from phenotype file"),
  make_option(c("--ld","-L"),type="character", help="Plink binary file prefix for the reference file used for LD calculation. 
                    If not provided, will use the target genotype for the LD calculation. Can also use multiple chromosome
                    plink file. Please see --target for more information."),
  make_option(c("--covar_header","-c"), type="character", help="Header of covariates. If not provided, will use all variable in the 
                    covariate file as the covarite."),
  make_option(c("--covar_file","-C"), type="character", help="Covarite file. Formate should be: ID Cov1 Cov2 ..."),
  make_option("--full",action="store_true", help="Also include the full model in the PRSice output"),
  make_option("--all", action="store_true", help="Output PRS for ALL threshold. Can only be used together with fastscore to avoid huge output files."),
  make_option("--no_regress", action="store_true", "Do not perform the regression analysis and simply output all PRS. Can only be used together 
                    with fastscore to avoid huge output files. If you must, you can modify bar_levels to obtain the fine scale PRS outputs"),
  make_option(c("--out","-o"),type="character", help="Prefix of all output.", default="PRSice"),
  make_option(c("--lower", "-l"), type="numeric", help="The starting p-value threshold. ", default=0.0001),
  make_option(c("--upper", "-u"), type="numeric", help="The final p-value threshold.", default=0.5),
  make_option(c("--interval", "-i"), type="numeric", help="The step size of the threshold.", default=0.00005),
  make_option("--fastscore", action="store_true", help="Calculate the minimum amount of threshold as required by the bar_level option"),
  make_option("--chr", type="character", help="Column header of Chromosome <Required>"),
  make_option("--A1", type="character", help="Column header of Reference Allele <Required>"),
  make_option("--A2", type="character", help="Column header of Alternaative Allele"),
  make_option("--stat", type="character", help="Column header of test statistic <Required>"),
  make_option("--snp", type="character", help="Column header of SNP id"),
  make_option("--bp", type="character", help="Column header of SNP location"),
  make_option("--se", type="character", help="Column header of Standard Error"),
  make_option(c("--pvalue","-p"), help="Column header of p-value <Required> "),
  make_option("--index",action="store_true", "Indicate all the above options are providing the INDEX of the corresponding column. 
                    (Index should be 0-based). Useful when your base file each have a different header but the column 
                    index remains the same"),
  make_option("--clump_p", type="numeric", help="The p-value threshold use for clumping.", default=1),
  make_option("--clump_r2", type="numeric", help="The R2 threshold for clumping. Please note that as we did not implement the
                    maximum likelihood R2 calculation, the clumping result can differ slightly from plink.", default=0.1),
  make_option("--clump_kb", type="numeric", help="The distance for clumping in kb.", default=250),
  make_option(c("--bed", "-B"), type="character", help="Bed file containing the selected regions. Name of bed file will be used as the region identifier."),
  make_option(c("--gtf", "-g"), type="character", help="GTF file containing gene boundaries. Required when --msigdb is set."),
  make_option(c("--msigdb", "-m"), type="character", help="MSIGDB file containing the pathway information require the gtf file."),
  make_option("--gen_bed", action="store_true", help="Generate bed file of gene regions from the gtf file."),
  make_option("--proxy", type="numeric", help="Proxy threshold for index SNP to be considered as part of the region represented by the clumped SNPs.
                    e.g. --proxy 0.8 means the index SNP will represent the region of any clumped SNPs that has a R2 >= 0.8 with it even if 
                    it is not physically within these regions"),
  make_option("--prslice", type="numeric", help="Perform PRSlice where the whole genome is first cut into bin size specified by this option. PRSice 
                    will then be performed on each bin. Bins are then sorted according to the their R2. PRSice is then performed again 
                    to find the best bin combination. This cannot be performed together with PRSet"),
  make_option("--bar_levels", type="character", help="Level of barchart to be plotted. When fastscore is set, PRSice will 
                    only calculate the PRS for threshold within the bar level"),
  make_option(c("--thread","-T"), type="numeric", help="Number of thread use", default=1),
  make_option("--c_help", action="store_true", help="Print the help message from the c++ program instead"),
  make_option("--plot", action="store_true", help="Indicate only plotting is required"),
  make_option("--intermediate", type="character", help="Pefix of the intermediate files for plotting (e.g. ignore .prsice and .best). If not provided, will deduce the file prefix from --base and --target"),
  make_option(c("--quantile", "-q"), type="numeric", help="Number of quantiles to plot. 0 = Not producing the quantile plot", default=0),
  make_option(c("--quant_extract", "-e"), type="character", help="File contain sample id to be plot on a separated quantile e.g. extra quantile containing only these samples" ),
  make_option("--bar_level", type="character", help="barchar level used for plotting", default="0.001,0.05,0.1,0.2,0.3,0.4,0.5"),
  make_option("--quant_ref", type="numeric", help="Reference quantile for quantile plot"),
  make_option("--scatter_r2",action="store_true", help="y-axis of the high resolution scatter plot should be R2"),
  make_option("--bar_col_r2", action="store_true", help="Change the colour of bar to R2 instead of p-value"),
  make_option("--bar_col_low", type="character", help="Colour of the poorest predicting thresholds", default="dodgerblue"),
  make_option("--bar_col_high", type="character", help="Colour of the highest predicting thresholds", default="firebrick"),
  make_option("--prsice", type="character", help="Location of the PRSice binary"),
  make_option("--dir", type="character", help="Location to install ggplot. Only require if ggplot is not installed")
)

argv <- parse_args(OptionParser(option_list=option_list))
stop()

argv = commandArgs(trailingOnly = TRUE)
help=(sum(c("--help", "-h") %in%argv)>=1)
if(help){
  print.arg.parser(p)
  quit();
}
argv <- parse_args(p)
not_cpp <- c("help", "c_help", "plot", "quantile", "quant_extract", "intermediate","quant_ref", "scatter_R2","bar_col_r2","bar_col_low","bar_col_high", "prsice", "dir")




# CALL_PRSICE: Call the cpp PRSice if required
# To ensure the excutable is set correctly
if(!is.na(argv$prsice)){
  if(!startsWith(argv$prsice, "/") && !startsWith(argv$prsice, ".")){
    argv$prsice = paste("./", argv$prsice, sep="")
  }
}
if(argv$c_help){
  if(is.na(argv$prsice)){
    stop("Cannot use c_help without specifying the location of the PRSice binary!");
  }
  system(paste(argv$prsice," --help",sep=""))
  quit();
}

# We don't bother to check if the input is correct, the parameter should be checked by the c++ program
add_command <- function(input){
  if(length(input)==1){
    if(is.na(input)){
      return(NA);
    }else{
      return(input)
    }
  }else{
    return(paste(input,collapse=","))
  }
}
command=""
if(!argv$plot){
  for(i in names(argv)){
    # only need special processing for flags and specific inputs
    if(i=="index"){
      if(argv[[i]]) command = paste(command, " --",i,sep="")
    }else if(i=="gen_bed"){
      if(argv[[i]]) command = paste(command, " --",i,sep="")
    }else if(sum(i%in%not_cpp)!=0){
      # ignore 
    }else{
      temp = add_command(argv[[i]])
      if(!is.na(temp)){
        command = paste(command, " --",i, " ", temp, sep="")
      }
    }
  }
  if(nchar(command)==0){
    print.arg.parser(p)
    quit()
  }
  ret<-system2(paste(dir,"bin/PRSice",sep=""), command, stdout=TRUE, stderr=TRUE)
  if(attr(ret, "status")!=0){
    quit();
  }
}




# PLOTTING: Here contains all the function for plotting
# quantile_plot: plotting the quantile plots
quantile_plot <- function(PRS, PRS.best, pheno, prefix, argv){
  extract = NULL
  if(!is.na(argv$quant_extract)){
    extract = fread(argv$quant_extract, header=F, data.table=F)
  }
  quants <- as.numeric(cut(PRS.best[,2], breaks = quantile(PRS.best[,2], probs = seq(0, 1, 1/argv$quantile)), include.lowest=T))
  num_quant <- argv$quantile
  if(!is.null(extract)){
    quants[PRS.best[,1]%in%extract$V2] <- num_quant+1
    num_quant<-num_quant+1;
  }
  quant.ref <- ceiling(argv$quantile/2)
  if(!is.na(argv$quant_ref)){
    quant.ref <- argv$quant_ref;
    if(quant.ref > argv$quantile){
      quant.ref <- ceiling(argv$quantile/2)
      cat(paste("WARNING: reference quantile", quant.ref, "is greater than number of quantiles", argv$quantile, "\n Using middle quantile by default"))
    }
  }
  
  quants <- factor(quants, levels = c(quant.ref, seq(1, num_quant, 1)[-quant.ref]))
  pheno$quantile <- quants
  if(ncol(pheno)>=3){
    pheno <- pheno[,c(colnames(pheno)[2],"quantile",colnames(pheno)[3:(ncol(pheno)-1)])]
  }else{
    pheno <- pheno[,c(colnames(pheno)[2],"quantile")]
  }
  family <- gaussian
  if(sum(unique(pheno[,1])%in%c(0,1)) ==2){
    # When only contain 0 and 1, we will consider it as binary
    # Someone will have to be very unlucky to have a quantitative trait as exactly 0 and 1 for all samples
    family <- binomial
  }
  reg <- summary(glm(Pheno ~ ., family, data = pheno))
  coef.quantiles <- reg$coefficients[1:num_quant,1]
  ci.quantiles.u <- reg$coefficients[1:num_quant,1] + (1.96*reg$coefficients[1:num_quant,2])
  ci.quantiles.l <- reg$coefficients[1:num_quant,1] - (1.96*reg$coefficients[1:num_quant,2])
  coef.quantiles[1] <- 0 
  ci.quantiles.u[1] <- 0
  ci.quantiles.l[1] <- 0
  quantiles.for.table <- c(quant.ref, seq(1, num_quant, 1)[-quant.ref])
  quantiles.df <- data.frame(coef.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
  names(quantiles.df) <- c("Coef", "CI.U", "CI.L", "DEC")
  quantiles.df$Group=0;
  if(!is.null(extract)){
    quantiles.df$Group[max(quantiles.df$DEC)]=1
  }
  quantiles.df$Group <- factor(quantiles.df$Group, levels=c(0,1))
  quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
  quantiles.plot <- ggplot(quantiles.df) + 
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5)) + 
    ylab("Change in Phenotype given score in quantiles") + 
    xlab("Quantiles for Polygenic Score") + 
    scale_x_continuous(breaks=seq(0, num_quant, 1)) +
    theme(axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"))
  if(is.null(extract)){
    quantiles.plot <- quantiles.plot+geom_point(aes(x = DEC, y = Coef), colour = "royalblue2", size=4) + 
      geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = Coef, x = DEC), colour = "royalblue2", size = 0.9)
  }else{
    quantiles.plot <- quantiles.plot+geom_point(aes(x = DEC, y = Coef,color=Group), size=4) + 
      geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = Coef, x = DEC, color=Group), size = 0.9)+
      scale_colour_manual(values=c("#0072B2", "#D55E00"))
  }
  ggsave(paste(prefix, "QUANTILES_PLOT.png", sep = "_")) 
}

high_res_plot <- function(PRS, prefix, argv){
  # we will always include the best threshold
  barchart.levels <- c(argv$bar_level, PRS$Threshold[which.max(PRS$R2)])
  barchart.levels <- sort(unique(barchart.levels),decreasing=F)
  # As the C++ program will skip thresholds, we need to artificially add the correct threshold information
  threshold_presented <- barchart.levels %in% PRS$Threshold
  for(i in 1:length(threshold_presented)){
    if(!threshold_presented[i]){
      barchart.levels[i] <- tail(PRS[PRS$Threshold<barchart.levels[i],],n=1)$Threshold
    }
  }
  # Need to also plot the barchart level stuff with green
  ggfig.points <- NULL
  if(argv$scatter_r2){
    ggfig.points <- ggplot(data=PRS, aes(x = Threshold, y = R2))+
      geom_line(aes(Threshold,  R2), colour = "green", data = PRS[with(PRS, Threshold %in% barchart.levels) , ] )+ 
      geom_hline(yintercept=max(PRS$R2),colour="red")+
      ylab(expression(paste("PRS model fit:  ", R^2, sep = " ")))
  }else{
    ggfig.points <- ggplot(data=PRS, aes(x = Threshold, y = -log10(P)))+
      geom_line(aes(Threshold,  -log10(P)), colour = "green", data = PRS[with(PRS, Threshold %in% barchart.levels) , ] )+ 
      geom_hline(yintercept=max(-log10(PRS$P)),colour="red")+
      ylab(bquote(PRS~model~fit:~italic(P)-value~(-log[10])))
  }
  ggfig.points <- ggfig.points + geom_point() + geom_line() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5), axis.line.x = element_line(color="black"),
                axis.line.y = element_line(color="black"))+
    xlab(expression(italic(P)-value~threshold~(italic(P)[T])));
   ggsave(paste(prefix,"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
}

bar_plot<-function(PRS, prefix, argv){
  barchart.levels <- c(argv$bar_level, PRS$Threshold[which.max(PRS$R2)])
  barchart.levels <- sort(unique(barchart.levels),decreasing=F)
  # As the C++ program will skip thresholds, we need to artificially add the correct threshold information
  threshold_presented <- barchart.levels %in% PRS$Threshold
  for(i in 1:length(threshold_presented)){
    if(!threshold_presented[i]){
      barchart.levels[i] <- tail(PRS[PRS$Threshold<barchart.levels[i],],n=1)$Threshold
    }
  }
  output <- PRS[PRS$Threshold %in% barchart.levels,]
  output$print.p[round(output$P, digits = 3) != 0] <- round(output$P[round(output$P, digits = 3) != 0], digits = 3)
  output$print.p[round(output$P, digits = 3) == 0] <- format(output$P[round(output$P, digits = 3) == 0], digits=2)
  output$print.p <- sub("e", "*x*10^", output$print.p)
  
  ggfig.plot <- ggplot(data=output)
  if(!argv$bar_col_r2){
    ggfig.plot <- ggfig.plot + geom_bar(aes(x = factor(Threshold), y = R2, fill = factor(Threshold)), stat="identity") +
      scale_fill_brewer(palette="YlOrRd", name = expression(italic(P)-value~threshold))
  }
  if(argv$bar_col_r2){
    ggfig.plot <- ggfig.plot + geom_bar(aes(x = factor(Threshold), y = R2, fill = -log10(P)), stat="identity") +     
      scale_fill_gradient(low= argv$bar_col_low, high= argv$bar_col_high, name =bquote(atop(-log[10]~model,italic(P)-value),))  
  }
  ggfig.plot <- ggfig.plot + geom_text(aes(x = factor(Threshold), y = R2, label = paste(print.p)), vjust = -1.5, hjust = 0, angle = 45, cex = 2.8, parse=T) +
    scale_y_continuous(limits = c(0, max(output$R2)*1.25)) +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5) , axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"))+
    xlab(expression(italic(P)-value~threshold~(italic(P)[T]))) + 
      ylab(expression(paste("PRS model fit:  ",R^2)))
    ggsave(paste(prefix, "_BARPLOT_", Sys.Date(), ".png", sep = ""))		
}  
# run_plot: The function used for calling different plotting functions
run_plot<-function(prefix, argv){
  PRS <- fread(paste(prefix,".prsice",sep=""), header=T,data.table=F)
  PRS.best <- fread(paste(prefix, ".best",sep=""), header=T,data.table=F)
  pheno <- fread(paste(prefix, ".pheno", sep=""), header=T, data.table=F)
  if(argv$quantile > 0){
    # Need to plot the quantile plot
    quantile_plot(PRS, PRS.best, pheno, prefix, argv)
  }
  # Now perform the barplotting
  if(!argv$fastscore){
    high_res_plot(PRS, prefix, argv)
  }
  bar_plot(PRS,prefix, argv)
}

# CALL PLOTTING FUNCTION: Process the input names and call the actual plotting function
if(!is.na(argv$intermediate)){
  # File name provided
  run_plot(argv$intermediate, argv)
}else{
  # we need to deduce the file name
  writeLines(strwrap("WARNING: Using this method, we will only perform plotting to the base region. If you are using PRSlice, please specify the name of the desired intermediate file",width=80))
  if(sum(is.na(argv$base))!=length(argv$base)){
    if(!is.na(argv$target)){
      print(argv$base)
      for(b in argv$base){
        print(argv$base)
        if(sum(!is.na(argv$pheno_col)) != 0){
          for(t in argv$pheno_col){
            # no need to check out, as it has default (PRSice)
            run_plot(paste(argv$out, b,t,"base",sep="."), argv);
          }
        }else{
          run_plot(paste(argv$out, b,"base",sep="."), argv);
        }
      }
    }else{
      message <- "NA in --target"
      if(argv$plot){
        message<-paste(message, "You'll need to either provide the intermediate prefix or all the target/base name for plot only")
      }
      stop(message);
    }
  }else{
    message <- "NA in --base"
    if(argv$plot){
      message<-paste(message, "You'll need to either provide the intermediate prefix or all the target/base name for plot only")
    }
    stop(message);
  }
}