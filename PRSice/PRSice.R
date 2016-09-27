# Here is the guide to this protentially long R code
# To go to each section, just search for the corresponding header as stated here
# The code structure are as follow
# INSTALL_PACKAGE
# - Contains functions responsible for installing all required packages
# COMMAND_PARSE
# - Functions and scripts to process all the command line arguments



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

dir.create(file.path("lib"), showWarnings = FALSE)
.libPaths(c(.libPaths(), "./lib"))
libraries <- c( "ggplot2", "argparser")
for(library in libraries)
{
  if(!UsePackage(library))
  {
    stop("Error: ", library, " cannot be load nor install!")
  }
}



# COMMAND_PARSE:  Functions and scripts to process all the command line arguments

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
    width=2*max_argument_length
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


p <- arg_parser_self("PRSice: Polygenic Risk Score software")
p <- add_argument(p, "--base", short="-b", nargs=Inf, help="Base association files, can input multiple times")
p <- add_argument(p, "--target", short="-t", nargs=Inf, help="Plink binary file prefix for target files. Can input multiple times")
p <- add_argument(p, "--binary_target", nargs=Inf, help="Indication of whether if the target is binary or not.Should be of the same  length as target")
p <- add_argument(p, "--superLongFlagHahaha", flag=T, help="Indication of whether if the target is binary or not.Should be of the same  length as target")
p <- add_argument(p, "positionFlagHahahajusttestingagain", help="Indication of whether if the target is binary or not.Should be of the same  length as target")

print.arg.parser(p)