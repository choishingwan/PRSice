# Here is the guide to this protentially long R code
# To go to each section, just search for the corresponding header as stated here
# The code structure are as follow
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
if(Sys.info()[1]=="Windows"){
  print("Window not supported, because of slash...")  
  print("They use backward slash, which will cause problems in the script as it is an escape character")
}

initial.options <- commandArgs(trailingOnly = FALSE)
if(length(commandArgs(trailingOnly = TRUE))==0){
  stop("Please use --help for the help messages")
}

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
dir=strsplit(script.name, "/")
dir = paste(dir[[1]][-length(dir[[1]])], "/", sep="",collapse="/")

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

UsePackage <- function(package)
{
  if(!InstalledPackage(package))
  {
    suppressMessages(suppressWarnings(install.packages(package, lib="./lib",repos="http://cran.rstudio.com/")))
    if(!InstalledPackage(package)) return(FALSE)
  }
  return(TRUE)
}

dir.create(file.path(dir,"lib"), showWarnings = FALSE)
.libPaths(c(.libPaths(), paste(dir,"lib",sep="")))
libraries <- c( "ggplot2", "argparser", "data.table")
for(library in libraries)
{
  if(!UsePackage(library))
  {
    stop("Error: ", library, " cannot be load nor install!")
  }
}



# COMMAND_FUNC:  Functions required for command line argument parsing

# Self defined arg_parser object. Difference from the default = removed place holder and the -x flag
arg_parser_self <- function(description, name = NULL)
{
  if (is.null(name)) {
    prefix <- "--file="
    name <- sub(prefix, "", grep(paste(prefix, "(.+)", sep = ""), commandArgs(), value = TRUE))
  }
  if (length(name) == 0) name <- "<script>"
  parser <- structure(list(name = name, description = description), class = "arg.parser")
  parser <- add_argument(parser, "--help", "show this help message and exit", 
                         flag = TRUE)
  parser    
}

# Self defined print function, should look much better than the default one
print.arg.parser<- function (x, width=NULL,...) 
{
  parser <- x
  opt.args <- parser$args[parser$is.opt.arg]
  
  max_name_length <- max(nchar(paste(parser$args[!is.na(parser$shorts)],parser$shorts, sep=", ")), nchar(parser$args[is.na(parser$shorts)]))+1
  max_argument_length <- max_name_length+max(nchar(toupper(sub("^-+", "", parser$args[parser$is.opt.arg]))))+1
  if(is.null(width)){
    width=max(80, 2*max_argument_length)
  }
  usage<-c("usage: ", parser$name, 
               paste(sub("^(.*)$", "[\\1]", parser$args[parser$is.flag])),
               paste(sub("^(.*)$", "[\\1 ", opt.args), toupper(sub("^--(.*)$", "\\1]", opt.args)), sep = ""),
               paste(parser$args[parser$is.req.arg]))
  current_message = paste(usage[1:2],sep=" ",collapse = "")
  pad = paste(rep(" ",nchar(current_message)),collapse = "");
  for(i_usage in 3:length(usage)){
    if(nchar(paste(current_message, usage[i_usage],sep=" ",collapse = ""))> width){
      writeLines(current_message)
      current_message = paste(pad,usage[i_usage],sep=" ",collapse = "");
    }else{
      current_message=paste(current_message, usage[i_usage], sep=" ",collapse = "")
    }
  }
  if(current_message!=pad){
    writeLines(current_message)
  }
  writeLines(strwrap(parser$description,width=width))
  
  
  if (sum(parser$is.req.arg) > 0) {
    message("positional arguments:")
    for (i_reg_arg in which(parser$is.req.arg)) {
      parser.frag = strwrap(parser$helps[i_reg_arg], width-6-max_argument_length)
      writeLines(paste("  ", parser$args[i_reg_arg], paste(rep(" ", max_argument_length-nchar(parser$args[i_reg_arg])), collapse=""), "    ",parser.frag[1],sep="",collapse=""))
      if(length(parser.frag)>1){
        for(i_frag in 2:length(parser.frag)){
          writeLines(paste(paste(rep(" ", 6+max_argument_length),collapse = ""), parser.frag[i_frag], collapse="",sep=""))
        }
      }
    }
  }
  
  message("")
  if (sum(parser$is.flag) > 0) {
    message("flags:")
    for (i_flag in which(parser$is.flag)) {
      if (parser$args[i_flag] == "--") 
        next
      if (is.na(parser$shorts[i_flag])) {
        arg.name <- parser$args[i_flag]
      }
      else {
        arg.name <- paste(parser$shorts[i_flag], parser$args[i_flag],  sep = ", ")
      }
      if(nchar(arg.name)<max_argument_length){
        arg.name <- paste(arg.name, paste(rep(" ", max_argument_length-nchar(arg.name)), collapse = ""), sep="")
      }
      arg.help <- argparser:::make_arg_help(parser, i_flag)
      arg.help.frag = strwrap(arg.help, width-6-max_argument_length)
      writeLines(paste("  ",arg.name, "    ", arg.help.frag[1],collapse = "",sep=""))
      if(length(arg.help.frag)>1){
        for(i_flag in 2:length(arg.help.frag)){
          writeLines(paste(paste(rep(" ", 6+max_argument_length),collapse = ""), arg.help.frag[i_flag], collapse="",sep=""))
        }
      }
    }
  }
  
  message("")
  if (sum(parser$is.opt.arg) > 0) {
    message("optional arguments:")
    for (i_opt_arg in which(parser$is.opt.arg)) {
      if (is.na(parser$shorts[i_opt_arg])) {
        arg.name <- parser$args[i_opt_arg]
      }
      else {
        arg.name <- paste(parser$shorts[i_opt_arg], parser$args[i_opt_arg],  sep = ", ")
      }
      
      if(nchar(arg.name)<max_name_length){
        arg.name <- paste(arg.name, paste(rep(" ", max_name_length-nchar(arg.name)), collapse = ""), sep="")
      }
      arg.name <- paste(arg.name, toupper(sub("^-+", "", parser$args[i_opt_arg])))
      if(nchar(arg.name)<max_argument_length){
        arg.name <- paste(arg.name, paste(rep(" ", max_argument_length-nchar(arg.name)), collapse = ""), sep="")
      }
      arg.help <- argparser:::make_arg_help(parser, i_opt_arg)
      arg.help.frag = strwrap(arg.help, width-6-max_argument_length)
      writeLines(paste("  ",arg.name, "    ", arg.help.frag[1],collapse = "",sep=""))
      if(length(arg.help.frag)>1){
        for(i_flag in 2:length(arg.help.frag)){
          writeLines(paste(paste(rep(" ", 6+max_argument_length),collapse = ""), arg.help.frag[i_flag], collapse="",sep=""))
        }
      }
    }
  }
}

# COMMAD_BUILD: Building the command line parser using argparser
# TODO: Write our own package to better handle these parameters. Currently they are all over the places
p <- arg_parser_self("PRSice: Polygenic Risk Score software")
p <- add_argument(p, "--base", short="-b", nargs=Inf, help="Base association files. User can provide multiple base files.")
p <- add_argument(p, "--target", short="-t", nargs=Inf, help="Plink binary file prefix for target files. User can provide multiple target files. Currently only support plink binary input. Does not support multi-chromosome input")
p <- add_argument(p, "--binary_target", nargs=Inf, help="Indication of whether binary target is provided. Should be of the same length as target")
p <- add_argument(p, "--beta", nargs=Inf, help="Indication of whether the test statistic is beta instead of OR. Should be of the same length as base")
p <- add_argument(p, "--pheno_file", short="-f", nargs=Inf, help="Phenotype file(s) containing the target phenotypes. If provided, the fam file of the target is ignored. This should be the same line as target (If you want to use phenotype file, you must use it for ALL target")
p <- add_argument(p, "--ld", short="-L", help="Plink binary file prefix for the reference file used for LD calculation. If not provided, will use the target genotype for the LD calculation")
p <- add_argument(p, "--covar_header", short="-c", help="Header of covariates, if not provided, will use all variable in the covariate file as the covarite. Should be comma separated")
p <- add_argument(p, "--covar_file", short="-C", help="Covariate file. Format should be: ID Cov1 Cov2 Must contain a header");
p <- add_argument(p, "--ancestry", short="-a", help="NOT DEVELOPED YET");
p <- add_argument(p, "--out", short="-o", help="The prefix of all output", default="PRSice");
p <- add_argument(p, "--lower", short="-l", help="The starting p-value threshold", default=0.0001);
p <- add_argument(p, "--upper", short="-u", help="The final p-value threshold", default=0.5);
p <- add_argument(p, "--interval", short="i", help="The step size of the threshold", default=0.00005);
p <- add_argument(p, "--chr", help="Column header of Chromosome", default="CHR");
p <- add_argument(p, "--A1", help="Column header of Reference Allele", default="A1");
p <- add_argument(p, "--A2", help="Column header of Alternaative Allele", default="A2");
p <- add_argument(p, "--stat", help="Column header of test statistic, either BETA or OR", default="OR");
p <- add_argument(p, "--snp", help="Column header of SNP id", default="SNP");
p <- add_argument(p, "--bp", help="Column header of SNP location", default="BP");
p <- add_argument(p, "--se",help="Column header of Standard Error", default="SE");
p <- add_argument(p , "--pvalue", short="-p", help="Column head of p-value", default="P");
p <- add_argument(p, "--index", flag=T, help="If the base file doesn't contain a header, you can use this option, which essentially state that all the provided \"headers\" are INDEX of the corresponding column. (Index should be 0-based)");
p <- add_argument(p, "--clump-p", help="The p-value threshold use for clumping", default=1);
p <- add_argument(p, "--clump_r2", help="The R2 threshold for clumping. Please note that as we did not implement the maximum likelihood R2 calculation, the clumping result can differ slightly from plink.", default=0.1);
p <- add_argument(p, "--clump_kb", help="The distance for clumping in kb.", default=250);
p <- add_argument(p, "--bed", short="-B", nargs=Inf, help="Bed file containing the selected regions. Name of bed file will be used as the region identifier");
p <- add_argument(p, "--gtf", short="-g", help="GTF file containing gene boundaries. Required when --msigdb is set");
p <- add_argument(p, "--msigdb", short="-m", help="MSIGDB file containing the pathway information require the gtf file");
p <- add_argument(p, "--gen_bed", flag=T, help="Generate bed file of gene regions from the gtf file");
p <- add_argument(p, "--proxy", help="Proxy threshold for index SNP to be considered as part of the region represented by the clumped SNPs. e.g. --proxy 0.8 means the index SNP will represent the region of any clumped SNPs that has a R2 >= 0.8 with it");
p <- add_argument(p, "--thread", short="-T", help="Number of thread used", default=1);
p <- add_argument(p, "--c_help", flag=T, help="Print the help message from the c++ program instead");
p <- add_argument(p, "--plot", flag=T, help="Indicate whether only plotting is required");
p <- add_argument(p, "--intermediate", help="Pefix of the intermediate files for plotting (e.g. ignore .prsice and .best). If not provided, will deduce the file prefix from --base and --target")
p <- add_argument(p, "--quantile", short="-q", help="Number of quantiles to plot. 0 = Not producing the quantile plot", default=0);
p <- add_argument(p, "--quant_extract", short="-e", help="File contain sample id to be plot on a separated quantile e.g. extra quantile containing only these samples" )
p <- add_argument(p, "--bar_level", help="barchar level used for plotting", default=c(0.001,0.05,0.1,0.2,0.3,0.4,0.5))
p <- add_argument(p, "--quant_ref", help="Reference quantile for quantile plot")
p <- add_argument(p, "--scatter_r2", flag=T, help="y-axis of the high resolution scatter plot should be R2")
p <- add_argument(p, "--fastscore", flag=T, help="Calculate the minimum amount of threshold as indicated by the bar_level option")
p <- add_argument(p, "--bar_col_r2", flag=T, help="Change the colour of bar to R2 instead of p-value")
p <- add_argument(p, "--bar_col_low", help="Colour of the poorest predicting thresholds", default="dodgerblue")
p <- add_argument(p, "--bar_col_high", help="Colour of the highest predicting thresholds", default="firebrick")

argv = commandArgs(trailingOnly = TRUE)
help=(sum(c("--help", "-h") %in%argv)>=1)
if(help){
  print.arg.parser(p)
  quit();
}
argv <- parse_args(p)

not_cpp <- c("help", "c_help", "plot", "quantile", "quant_extract", "intermediate","quant_ref", "scatter_R2","bar_col_r2","bar_col_low","bar_col_high")
# CALL_PRSICE: Call the cpp PRSice if required
if(argv$c_help){
  system(paste(dir,"bin/PRSice --help",sep=""))
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
    if(sum(is.na(argv$target))!=length(argv$target)){
      for(b in argv$base){
        for(t in argv$target){
          # no need to check out, as it has default (PRSice)
          run_plot(paste(argv$out, b,t,"base",sep="."), argv);
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