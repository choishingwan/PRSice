library(batch)
start.time <- proc.time()[3]
options(echo = FALSE)
options(warn=-1)

cat(" ################################# \n # \n # \n # \n # \n # PRSice: Polygenic Risk Score software \n # \n # Jack Euesden, Cathryn M. Lewis, Paul F. O'Reilly 2014 \n # \n # \n # If you use PRSice in published work, please cite: \n # \n # \"PRSice: Polygenic Risk Score software\" \n # Euesden, Lewis, O'Reilly, Bioinformatics (2015) 31 (9):1466-1468 \n # \n # \n # \n #  \n ################################# \n")



#################################
#
#  Default options
#
#################################
								
## essential
target <-  NA
base <-   NA

## preferable
plink <-  NA
order.cols <- "SNP,CHR,BP,A1,A2,OR,SE,P"
supplied.order <- F

# phenotype options
pheno.file <-   NA
binary.target <-  T

# covariate options
covary <- T
covariates <- "C1,C2"
user.covariate.file <-  NA
ancestry.dim <- "MDS"

# graphical parameters
ggfig <-  T
barchart.levels <- "0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5"
barpalatte <- "YlOrRd"
best.thresh.on.bar <- F
scatter.R2 <- F						
figname <- "PRSice"
bar.col.is.pval <- T
bar.col.is.pval.lowcol <- "dodgerblue"
bar.col.is.pval.highcol <- "firebrick"

# clumping							
clump.snps <- T
clump.p1 <- 1
clump.p2 <- 1
clump.r2 <- 0.1
clump.kb <- 250

# pruning
prune.snps <- F
prune.kb.wind <- 50
prune.kb.step <- 2
prune.kb.r2 <- 0.8

							
# high density scoring							
slower <- 0.0001
supper <- 0.5
sinc <-  0.00005
fastscore <- F
							
# dosage
dosage <- F
dosage.format <- "gen"
dos.skip0 <- 1
dos.skip1 <- 1
dos.coding <- 1
dos.format <- 3
dos.sep.fam <- NA
dos.fam.is.samp <- F
dos.impute2 <- F
dos.path.to.lists <- NA
dos.list.file <- NA

# alternate genotype formats
geno.is.ped <- F
geno.as.list <- F

# miscellaneous options
wd <-  "./"							
print.time <-  T
cleanup <- T
plinkpath <-  "./"
remove.mhc <- F
for.meta <- F
plink.silent <- T
allow.no.sex <- F
report.snps.best.model <- F

# A Better Coefficient of Determination
calculate.abc.r2 <- F
pi0 <- 0
n.ca.base <- NA
n.co.base <- NA
n.ca.targ <- NA
n.co.targ <- NA
prev.base <- NA
prev.targ <- NA

# quantiles
quantiles <- F
num.quantiles <- 5
quant.ref <- NA

mend.score <- F
mend.score.len <- 100
score.at.1 <- F
report.individual.scores <- T
report.best.score.only <- T
no.regression <- F
debug.mode <- F

# sumsum
sumsum <- F
clump.ref <- NA
size.targ <- NA

# multiple phenotype options
multiple.target.phenotypes <- F
target.phenotypes <-  NA #  "V2,V3"
target.phenotypes.binary <- NA # "T,T"  ## NB:: 'QT' or 'BIN' 
multiple.base.phenotypes <- F
base.phenotypes.names <- NA  ## sub in for PHEN.NAME
heat.r2 <- F

# empirical p-value options:
n.emp.perms <- 10000
emp.pval <- F
report.perm.phen <- F
report.perm.pvals <- F
emp.alpha <- F

if(!fastscore){
  best.thresh.on.bar <- T
}

if(dosage){
  quantiles <- F
}



#cat(" ################################# \n # \n #  Begin Script \n # \n ################################# \n")


if(Sys.info()[1] == "Darwin"){
	os <- "mac"
}
if(Sys.info()[1] == "Linux"){
	os <- "linux"
}
if(Sys.info()[1] == "Windows"){
	os <- "windows"
}

if(os == "windows"){
	print("ERROR: Windows not supported")
	quit()
}
													

if(dosage){
	sinc <- 0.001
	slower <- 0.001
	supper <- 0.5
}


cat(" ################################# \n # \n #  Read in Command Line Arguments & interpret \n # \n ################################# \n")


parseCommandArgs(evaluate=T)

supper <- supper+sinc
slower <- slower - sinc
covariates <- strsplit(covariates, split=",")[[1]]
order.cols <- strsplit(order.cols, split=",")[[1]]
barchart.levels <- as.numeric(strsplit(barchart.levels, split=",")[[1]])
ggfig <- as.logical(ggfig)
covary <- as.logical(covary)
fastscore <- as.logical(fastscore)
print.time <- as.logical(print.time)
supplied.order <- as.logical(supplied.order)
binary.target <- as.logical(binary.target)
best.thresh.on.bar <- as.logical(best.thresh.on.bar)
cleanup <- as.logical(cleanup)
clump.snps <- as.logical(clump.snps)
prune.snps <- as.logical(prune.snps)
scatter.R2 <- as.logical(scatter.R2)
dosage <- as.logical(dosage)
dos.fam.is.samp <- as.logical(dos.fam.is.samp)
dos.impute2 <- as.logical(dos.impute2)
remove.mhc <- as.logical(remove.mhc)
bar.col.is.pval <- as.logical(bar.col.is.pval)
for.meta <- as.logical(for.meta)
geno.is.ped <- as.logical(geno.is.ped)
geno.as.list <- as.logical(geno.as.list)
mend.score <- as.logical(mend.score)
score.at.1 <- as.logical(score.at.1)
report.individual.scores <- as.logical(report.individual.scores)
report.best.score.only <- as.logical(report.best.score.only)
plink.silent <- as.logical(plink.silent)
no.regression <- as.logical(no.regression)
multiple.target.phenotypes <- as.logical(multiple.target.phenotypes)
multiple.base.phenotypes <- as.logical(multiple.base.phenotypes)
debug.mode <- as.logical(debug.mode)
quantiles <- as.logical(quantiles)
sumsum <- as.logical(sumsum)
calculate.abc.r2 <- as.logical(calculate.abc.r2)
heat.r2 <- as.logical(heat.r2)
emp.pval <- as.logical(emp.pval)
report.perm.phen <- as.logical(report.perm.phen)
report.perm.pvals <- as.logical(report.perm.pvals)
emp.alpha <- as.logical(emp.alpha)
report.snps.best.model <- as.logical(report.snps.best.model)

barchart.levels.old <- barchart.levels

if(multiple.base.phenotypes){
  base.phenotypes.names <- strsplit(gsub(" ", "", base.phenotypes.names), split=",")[[1]]
  if(is.na(base.phenotypes.names)){
    cat("ERROR: Please select base phenotypes to use, using base.phenotypes.names \n or set multiple.base.phenotypes F. \n Quitting")
    quit()
  }
}

if(multiple.target.phenotypes){
  target.phenotypes <- strsplit(gsub(" ", "", target.phenotypes), split=",")[[1]]  
  target.phenotypes.binary  <- as.logical(strsplit(gsub(" ", "", target.phenotypes.binary), split=",")[[1]])
  if(is.na(target.phenotypes)){
    cat("ERROR: Please select target phenotypes to use, using target.phenotypes \n or set multiple.target.phenotypes F. \n Quitting")
    quit()
  }
  if(is.na(target.phenotypes.binary)){
    cat("ERROR: Please specify whether target phenotypes are binary or QT \n or set multiple.target.phenotypes F. \n Quitting")
    quit()
  }
  if(length(target.phenotypes) != length(target.phenotypes.binary)){
  	cat("ERROR: Different number of target phenotype names and target phenotype types specified \n Check these lists are the same length \n Quitting")
  	quit()
  } 
}

if(!sumsum & quantiles){
 if(is.na(quant.ref)){
 	quant.ref <- ceiling(num.quantiles/2)
 }
 if(!is.na(quant.ref)){
 	if(quant.ref > num.quantiles){
   	  quant.ref <- ceiling(num.quantiles/2)
 	  cat(paste("WARNING: reference quantile", quant.ref, "is greater than number of quantiles", num.quantiles, "\n Using middle quantile by default"))
    }
  }
}

if(calculate.abc.r2 & binary.target){
  n1 <- n.ca.base + n.co.base
  n2 <- n.ca.targ + n.co.targ 
  sampling1 <- n.ca.base / n1
  sampling2 <- n.ca.targ / n2
  if(is.na(n.ca.base) | is.na(n.co.base) | is.na(n.ca.targ) | is.na(n.co.targ) ){
  	cat("ERROR: Need to supply information on \n Base and Target sizes \n to use R2 on Liability Scale \n Defaulting to Nagelkerke's Pseudo R2 \n")
    calculate.abc.r2 <- F
  }
  prevalence1 <- prev.base
  prevalence2 <- prev.targ
  if(is.na(prev.targ) | is.na(prev.base) ){
  	cat("ERROR: Need to supply information on \n Base and Target phenotype prevalences \n to use R2 on Liability Scale \n Defaulting to Nagelkerke's Pseudo R2 \n")
    calculate.abc.r2 <- F
  }
}

if(debug.mode){
  options(warn=0)
  plink.silent <- F
}

if(!is.na(pheno.file)){
  ext.phen <- T
}
if(is.na(pheno.file)){
  ext.phen <- F
}
if(!is.na(user.covariate.file)){
  covary <- T
}

if(!sumsum){
if(prune.snps & clump.snps){
  cat("ERROR:: Please select either clumping or pruning (or neither) \n Defaulting to Clumping \n")
  prune.snps <- F
}
}

if(multiple.target.phenotypes){
  ext.phen <- T
  cat("Using data on multiple phenotypes for target data \n")
}



######################################
#
#
#  Compile Functions from polygenescore.R (Dudbridge 2013): Begin
#
#
#######################################

## PolygeneScore from Dudbridge


polygenescore=function(n1,nsnp,vg1=1,n2=n1,vg2=vg1,corr=1,plower=0,pupper=1,weighted=T,alpha=0.05,binary=F,prevalence1=0.1,prevalence2=prevalence1,sampling1=prevalence1,sampling2=prevalence2,lambdaS1=NA,lambdaS2=NA,nullfraction=0,shrinkage=F,logrisk=F) {
    
    # variance of the phenotype
    varY1 = 1
    varY2 = 1
    if (binary) {
      varY1 = sampling1*(1-sampling1)
      varY2 = sampling2*(1-sampling2)
    }
    
    # sampling variance of each beta
    samplingVar = varY1/n1
    
    # conversion from lambdaS to vg
    if (!is.na(lambdaS1)) {
      if (logrisk) {
        vg1=2*log(lambdaS1)
      }
      else {
        t=-qnorm(prevalence1)
        t1=qnorm(1-lambdaS1*prevalence1)
        vg1 = 2*(t-t1*sqrt(1-(t^2-t1^2)*(1-t*prevalence1/dnorm(t))))/(dnorm(t)/prevalence1+t1^2*(dnorm(t)/prevalence1-t))
        if (vg1>1 | is.na(vg1)) vg1=1
      }
    }
    if (!is.na(lambdaS2)) {
      if (logrisk) {
        vg2=2*log(lambdaS2)
      }
      else {
        t=-qnorm(prevalence2)
        t1=qnorm(1-lambdaS2*prevalence2)
        vg2 = 2*(t-t1*sqrt(1-(t^2-t1^2)*(1-t*prevalence2/dnorm(t))))/(dnorm(t)/prevalence2+t1^2*(dnorm(t)/prevalence2-t))
        if (vg2>1 | is.na(vg2)) vg2=1
      }
    }
# variance of true betas
    betaVar = vg1/(nsnp*(1-nullfraction))
    betaVar2 = vg2/(nsnp*(1-nullfraction))
    
    # transform from liability scale to observed scale
    if (logrisk) {
      liab2obs1=prevalence1*sampling1*(1-sampling1)/prevalence1/(1-prevalence1)
      liab2obs2=prevalence2*sampling2*(1-sampling2)/prevalence2/(1-prevalence2)
    }
    else {
      liab2obs1=dnorm(qnorm(prevalence1))*sampling1*(1-sampling1)/prevalence1/(1-prevalence1)
      liab2obs2=dnorm(qnorm(prevalence2))*sampling2*(1-sampling2)/prevalence2/(1-prevalence2)
    }
    if (binary) betaVar = betaVar*liab2obs1^2
    if (binary) betaVar2 = betaVar2*liab2obs2^2
    
    shrink=1
    if (shrinkage) {
      shrink = 1-samplingVar/(betaVar*(1-nullfraction)+samplingVar)
    #  betaVar = betaVar*shrink^2
    #  betaVar2 = betaVar2*shrink^2
    #  samplingVar = samplingVar*shrink^2
    }
    
    # threshold on betahat based on its p-value
    betaHatThreshLo = -qnorm(plower/2)*sqrt(samplingVar)
    betaHatThreshHi = -qnorm(pupper/2)*sqrt(samplingVar)
    
    # expected number of selected SNPs
    betaHatSD = sqrt(betaVar+samplingVar)
    probTruncBeta = 2*nsnp*(1-nullfraction)*abs(pnorm(-betaHatThreshHi,sd=betaHatSD)-
                                                pnorm(-betaHatThreshLo,sd=betaHatSD))
    nullHatSD = sqrt(samplingVar)
    probTruncNull = 2*nsnp*nullfraction*abs(pnorm(-betaHatThreshHi,sd=nullHatSD)-
                                            pnorm(-betaHatThreshLo,sd=nullHatSD))
    
    # variance of the estimated gene score
    if (weighted) {
      if (plower==0) term1=0 else term1=betaHatThreshLo/betaHatSD*dnorm(betaHatThreshLo/betaHatSD)
      if (pupper==0) term2=0 else term2=betaHatThreshHi/betaHatSD*dnorm(betaHatThreshHi/betaHatSD)
      varBetaHat = betaHatSD^2*(1+(term1-term2)/(pnorm(betaHatThreshHi/betaHatSD)-pnorm(betaHatThreshLo/betaHatSD)))
      if (plower==0) term1=0 else term1=betaHatThreshLo/nullHatSD*dnorm(betaHatThreshLo/nullHatSD)
      if (pupper==0) term2=0 else term2=betaHatThreshHi/nullHatSD*dnorm(betaHatThreshHi/nullHatSD)
      varNullHat = samplingVar*(1+(term1-term2)/(pnorm(betaHatThreshHi/nullHatSD)-pnorm(betaHatThreshLo/nullHatSD)))
      varGeneScoreHat = varBetaHat*probTruncBeta+varNullHat*probTruncNull
    }
    else {
      varGeneScoreHat = probTruncBeta+probTruncNull
    }

    # covariance between Y2 and estimated gene score
    if (weighted) {
    # coefficient in SNPs with effects
      scoreCovariance = corr*sqrt(betaVar*betaVar2)/(betaVar+samplingVar)
    # covariance in SNPs with effects
      scoreCovariance = scoreCovariance*varBetaHat*probTruncBeta
    }
    else {
      scoreCovariance = 2*sqrt(betaVar2/betaVar)*corr*(1-nullfraction)*nsnp*
        integrate(discordantSign,0,Inf,sqrt(betaVar),betaHatThreshLo,betaHatThreshHi,sqrt(samplingVar),abs.tol=1e-12)$value
    }
    
    # Coefficient of determination!
    R2 = scoreCovariance^2/varGeneScoreHat/varY2
    # Non-centrality parameter!
    NCP=n2*R2/(1-R2)
    # Power!
    power=pchisq(qchisq(1-alpha,1),1,lower=F,ncp=NCP)
    
    thresholdDensity = dnorm(qnorm(prevalence2))/prevalence2
    caseMean = thresholdDensity*R2*varY2/liab2obs2^2
    caseVariance = R2*varY2/liab2obs2^2*(1-caseMean*(thresholdDensity+qnorm(prevalence2)))
    thresholdDensity = dnorm(qnorm(prevalence2))/(1-prevalence2)
    controlMean = -thresholdDensity*R2*varY2/liab2obs2^2
    controlVariance = R2*varY2/liab2obs2^2*(1+controlMean*(thresholdDensity-qnorm(prevalence2)))
    
    # debugging
    #print(c(probTruncBeta,probTruncNull,varGeneScoreHat,scoreCovariance,caseMean,controlMean,caseVariance,controlVariance))
    #print(varGeneScoreHat)
    
    # area under ROC curve!
    if (binary) {
      if (logrisk) {
        AUC=pnorm(sqrt(R2*(1-prevalence2)^2/sampling2/(1-sampling2)/2))
      }
      else {
        AUC = pnorm((caseMean-controlMean)/sqrt(caseVariance+controlVariance))
      }
      MSE=NULL
    }
    else {
      AUC = NULL
      MSE = 1+shrink^2*varGeneScoreHat-2*shrink*scoreCovariance
    }
    
    # R2 on liability scale for binary traits
    if (binary) R2=R2/liab2obs2^2*sampling2*(1-sampling2)
    
    return(list(R2=R2,NCP=NCP,p=pchisq(NCP+1,1,lower=F),power=power,AUC=AUC,MSE=MSE))
}


estimateVg2FromP=function(p,n1,nsnp,vg1=0,n2=n1,corr=1,plower=0,pupper=1,weighted=T,binary=F,prevalence1=0.1,prevalence2=prevalence1,sampling1=prevalence1,sampling2=prevalence2,lambdaS1=NA,lambdaS2=NA,nullfraction=0,shrinkage=F,logrisk=F) {
    obj1=function(vg2) {
        if (vg1==0) vg1here=vg2
        else vg1here=vg1
        (sqrt(polygenescore(n1=n1,nsnp=nsnp,vg1=vg1here,n2=n2,vg2=vg2,corr=corr,plower=plower,pupper=pupper,weighted=weighted,binary=binary,prevalence1=prevalence1,prevalence2=prevalence2,sampling1=sampling1,sampling2=sampling2,lambdaS1=lambdaS1,lambdaS2=lambdaS2,nullfraction=nullfraction,shrinkage=shrinkage,logrisk=logrisk)$NCP)-ncp)^2
    }
    ncp=qnorm(p/2,lower=F)
    vg=optimise(obj1,c(0,1))$minimum
    ncp=qnorm(.025,mean=qnorm(p/2,lower=F))
    vgLo=optimise(obj1,c(0,1))$minimum
    ncp=qnorm(.975,mean=qnorm(p/2,lower=F))
    vgHi=optimise(obj1,c(0,1))$minimum
    list(vg=vg,vgLo=vgLo,vgHi=vgHi)
}

######################################
#
#
#  Compile Functions from polygenescore.R (Dudbridge 2013): End
#
#
#######################################


if(!sumsum){
    for(basePhen in 1:length(base.phenotypes.names)){
	    output <- as.vector(1)
	
    	if(multiple.base.phenotypes){ 
            if(basePhen == 1){
                base.temp <- base
        	    base <- gsub("PHEN.NAME", base.phenotypes.names[basePhen], base)
      	    } 
    	    if(basePhen > 1){
      	        base <- gsub("PHEN.NAME", base.phenotypes.names[basePhen], base.temp)
      	    } 
            cat(paste("Constructing Polygenic Scores for ", base.phenotypes.names[basePhen], " Risk \n"))
        }
	
    	if(dosage){
	
	        cat(" ################################# \n # \n #  Dosage data defaults \n # \n ################################# \n")
	        if(dos.impute2){
        	    dosage.format <- "gen"
        	    dos.skip0 <- 1
        	    dos.skip1 <- 1
        	    dos.coding <- 1
        	    dos.format <- 3	
	        }
	        if(!dos.fam.is.samp){
	            target.fam <- dos.sep.fam
	        }
	        if(dos.fam.is.samp){
	            system(paste("tail -n +3 ", dos.sep.fam, " | awk '{print $1,$2,0,0,$4,$7}' > ./fam_from_sample", sep = ""))
	            target.fam <- "./fam_from_sample"
	        }
	        if(covary & !is.na(user.covariate.file)){
        	    print("ERROR: Cannot calculate covariates automatically for dosage files");quit()
	        }
	        paste("NB:: DOSAGE OPTION SELECTED. ANALYSIS MAY BE MUCH SLOWER")
	        if(!is.na(dos.path.to.lists) & is.na(dos.list.file)){
	            dos.list.file <- as.vector(1)
	            for(i in 1:22){
	                dos.list.file[i] <- paste(gsub(pattern="CHRNUM", replacement=seq(1,22,1)[i], x = dos.path.to.lists))
	            }
	            write.table(dos.list.file,"target_files",col.names=F,row.names=F,quote=F)
	            gen.name <- "target_files"
	            dos.list.file <- NA
	        }
	        if(!is.na(dos.list.file)){
	            gen.name <- dos.list.file
	        }
	        if(!is.na(target)){
	            gen.name <- target
	        }	
	    }
	
	
	    cat(" ################################# \n # \n #  Check options match \n # \n ################################# \n")
	
	    use.beta <- F
	
    	if(ancestry.dim != "MDS" & ancestry.dim != "PCA"){
    	    print("ERROR: cannot calculate that format of ancestry informative covariates")
    	    quit()
    	}
    	if(ancestry.dim == "PCA" & dosage){
    	    print("ERROR: plink-1.07 does not support principal components, plink-1.9 does not support dosages")
    	    quit()
    	}
	
    	if(is.na(target) & is.na(dos.path.to.lists) & is.na(dos.list.file)){
    	    print("ERROR: Please Supply a TARGET DATA SET");quit()
    	}
    	if(is.na(base)){
    	    print("ERROR: Please Supply a BASE DATA SET");quit()
    	}

	
	    if(is.na(plink)){
	        cat("ERROR: Please supply a path to PLINK. \n For Download Links, see www.PRSice.info \n NB: PLINK-1.07 is required for dosage data, PLINK2 is required for genotype data \n Quitting \n"); quit()
	    }
	    if(wd == "./"){
	        print("Using current directory as working directory")
	    }
	
	    setwd(wd)
	    mhc <- ""
	    if(remove.mhc){
		    mhc <- " --exclude mhc.txt range "
	    }
	    write.table(c("6 26000000 33000000 mhc"), "mhc.txt", col.names = F, row.names = F, quote = F)
	    if(ggfig){
	        library(ggplot2)
	        library(plyr)
    	}
	    if(binary.target){
	        library(fmsb)
	    }
	    options("scipen"=100,"digits"=4)
	    cat(" ################################# \n # \n #  Check base input format \n # \n ################################# \n")
	    if(basePhen == 1 ){
  	        if(plink.silent){
	            plink <- paste(plink, " --silent ")
	        }
	        if(allow.no.sex){
	            plink <- paste(plink, " --allow-no-sex ")
	        }
	    }
	    # default situation - assume a conventional header line to base
    	system(paste("head ", base, " > head_disc"))
    	head.disc <- read.table("head_disc", head = T)
    	col1 <- which(colnames(head.disc) == "SNP"); if(length(col1) == 0){print("ERROR: No SNP column in base data");quit()}
    	col2 <- which(colnames(head.disc) == "CHR"); if(length(col2) == 0){print("WARNING: No CHR column in base data")}
    	col3 <- which(colnames(head.disc) == "BP"); if(length(col3) == 0){print("WARNING: No BP column in base data")}
    	col4 <- which(colnames(head.disc) == "A1"); if(length(col4) == 0){print("ERROR: No A1 column in base data");quit()}
    	col5 <- which(colnames(head.disc) == "A2"); if(length(col5) == 0){print("WARNING: No A2 column in base data")}
    	col6 <- which(colnames(head.disc) == "OR"); if(length(col6) == 0){
	        col6 <- which(colnames(head.disc) == "BETA");print("NB: READING BETA NOT ODDS RATIO");use.beta <- T
	        if(length(col6) == 0){
		        print("ERROR: No OR column in base data");quit()
	        }
	    }	
	    col7 <- which(colnames(head.disc) == "SE"); if(length(col7) == 0){print("WARNING: No SE column in base data")}
	    col8 <- which(colnames(head.disc) == "P"); if(length(col8) == 0){print("ERROR: No P-val column in base data");quit()}
	    ## check for missing data in columns, reorder to account for this
	    if(length(col2) != 0){
	        chr.head <- "CHR"
	        chr.awk <- paste(",$", col2, sep = "")
	    }
    	if(length(col2) == 0){
    	    chr.head <- " "
    	    chr.awk <- ""
    	}
    	if(length(col3) != 0){
    	    bp.head <- "BP"
    	    bp.awk <- paste(",$", col3, sep = "")
    	}
    	if(length(col3) == 0){
    	    bp.head <- " "
    	    bp.awk <- ""
    	}
    	if(length(col5) != 0){
    	    a2.head <- "A2"
    	    a2.awk <- paste(",$", col5, sep = "")
    	}
    	if(length(col5) == 0){
    	    a2.head <- " "
    	    a2.awk <- ""
    	}
    	if(length(col7) != 0){
    	    se.head <- "SE"
    	    se.awk <- paste(",$", col7, sep = "")
    	}
    	if(length(col7) == 0){
    	    se.head <- " "
    	    se.awk <- ""
    	}
    	if(!use.beta){
    	    system(paste("awk '{print $",col1,chr.awk,bp.awk,",$",col4,a2.awk,",$",col6,se.awk,",$",col8,"}' ", base, " > reordered_base", sep  =""))
    	    system(paste("echo SNP",chr.head,bp.head,"A1",a2.head,"OR",se.head," P > HEADER"))
    	}
    	if(use.beta){
    	    system(paste("awk '{print $",col1,chr.awk,bp.awk,",$",col4,a2.awk,",$",col6,se.awk,",$",col8,"}' ", base, " > reordered_base", sep  =""))
    	    system(paste("echo SNP",chr.head,bp.head,"A1",a2.head,"BETA",se.head," P > HEADER"))
    	}
	
	    head.disc <- as.vector(t(read.table("HEADER", head = F)))
	
	
    	if(geno.is.ped & !geno.as.list){
    	    cat(" ################################# \n # \n #   Reformat Target \n # \n ################################# \n")
    	    system(paste(plink," --noweb --file", target,"--make-bed --out", target))
    	}
	    if(geno.is.ped & geno.as.list){
	        cat(" ################################# \n # \n #   Reformat Target \n # \n ################################# \n")
	        for(i in 1:22){
	            file.name.chr <- paste(gsub(pattern="CHRNUM", replacement=i, x = target))
	            system(paste(plink," --noweb --file", file.name.chr,"--make-bed --out", file.name.chr))
	        }
	    }
	
	    if(!dosage & !geno.as.list){
	        system(paste("awk '{print $",which(head.disc=="SNP"),"}' reordered_base | sort -k1,1 > base_SNPS", sep=""))
	        if(as.numeric(gsub(" ", "", system(paste("awk '{print $2}' ",target, ".bim  | sort -k1,1 | join -1 1 -2 1 ``-'' base_SNPS | wc -l",sep = ""),intern=T)))==0){
	            cat("ERROR: No SNP's present in both base and target data set. Exiting. \n"); quit()
	        }
	    }
	    if(!dosage & geno.as.list){
	        system(paste("awk '{print $",which(head.disc=="SNP"),"}' reordered_base | sort -k1,1 > base_SNPS", sep=""))
	        if(as.numeric(gsub(" ", "", system(paste(gsub("CHRNUM", "*", paste("cat ", target, "*.bim | awk '{print $2}'  | sort -k1,1 | join -1 1 -2 1 ``-'' base_SNPS | wc -l",sep = ""))),intern=T)))==0){
	            cat("ERROR: No SNP's present in both base and target dataset. Exiting. \n"); quit()
	        }
	    }
	    if(!dosage){
	        cat(" ################################# \n # \n #   Remove Ambiguous SNPs \n # \n ################################# \n")
	        if(!use.beta){
	            system(paste("tail -n +2 reordered_base | awk '{print $",which(head.disc=="SNP"),
					",$",which(head.disc=="A1"),",log($",which(head.disc=="OR"),")}' | sort -k1,1 > temp.raw", sep = ""))
	        }
	        if(use.beta){
	            system(paste("tail -n +2 reordered_base | awk '{print $",which(head.disc=="SNP"),
					",$",which(head.disc=="A1"),",$",which(head.disc=="BETA"),"}' | sort -k1,1 > temp.raw", sep = ""))
	        }
	
	        if(!geno.as.list){ 
	            system(capture.output(cat("awk '{print $2,$5,$6}' ", target,".bim | sort -k1,1 | join -1 1 -2 1   temp.raw   ``-'' | awk '($4 == \"A\" && $5 == \"T\" ||$4 == \"T\" && $5 == \"A\" || $4 == \"C\" && $5 == \"G\"  || $4 == \"G\" && $5 == \"C\"  ){print $1}' > synonymous_snps  \n ", sep = ""))[1])
	            system(paste( plink, "        --noweb --bfile ", target," --exclude synonymous_snps  --make-bed   --out non_synonymous_snps_only", sep = ""))
	        }
	        if(geno.as.list){
	  	        for(chrnum in 1:22){
	  	            file.name.chr <- paste(gsub(pattern="CHRNUM", replacement=chrnum, x = target))
	                system( capture.output(cat("awk '{print $2,$5,$6}' ", file.name.chr,".bim | sort -k1,1 | join -1 1 -2 1   temp.raw   ``-'' | awk '($4 == \"A\" && $5 == \"T\" ||$4 == \"T\" && $5 == \"A\" || $4 == \"C\" && $5 == \"G\"  || $4 == \"G\" && $5 == \"C\"  ){print $1}' > synonymous_snps_",chrnum,"  \n ", sep = ""))[1])
	                system(paste( plink, "        --noweb --bfile ", file.name.chr," --exclude synonymous_snps_",chrnum,"  --make-bed   --out non_synonymous_snps_only_",chrnum, sep = ""))
	            }
	            system("cat non_synonymous_snps_only_*.bim  > non_synonymous_snps_only.bim")  
	        }
	    }
	    if(clump.snps){
	        cat(" ################################# \n # \n #   Clump \n # \n ################################# \n")
	    }
	
	    if(dosage){
	        system(paste("tail -n +2 reordered_base   > cleaned_base ", sep = ""))
	    }
	    if(!dosage){
	        system(paste("awk '{print $2}' non_synonymous_snps_only.bim | sort -k 1,1 >  ./TARGET_SNPs", sep = ""))
	        system(paste("tail -n +2 reordered_base | sort -k1,1  | join -1 1 -2 1 ``-'' ./TARGET_SNPs   > cleaned_base ", sep = ""))
	    }
	
	    system("cat HEADER cleaned_base > cleaned_base.assoc")
	
	    if(clump.snps & !dosage & !(geno.as.list)){
            ## clump and extract independant snps
			system(paste(plink," --noweb --bfile non_synonymous_snps_only --clump cleaned_base.assoc --clump-p1 ",clump.p1," --clump-p2 ",clump.p2," --clump-r2 ",clump.r2," --clump-kb ",clump.kb," --out cleaned_base", sep=" "))
			system("tail -n +2 cleaned_base.clumped | awk '{print $3}'  | sort -k1,1 | awk '($1 != \"\"){print}'  > LE_SNPs")
			## linkage equilibrium snps
			if(!use.beta){
			    system(paste("tail -n +2 cleaned_base.assoc | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' | awk '{print $",which(head.disc=="SNP"),
					",$",which(head.disc=="A1"),",log($",which(head.disc=="OR"),")}' > rawfile.raw", sep = ""))
				system(paste("awk '{print $",which(head.disc=="SNP"),",$",which(head.disc=="P"),"}' cleaned_base.assoc | tail -n +2 | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' > rangelist_ranges"))
			}
			if(use.beta){
				system(paste("tail -n +2 cleaned_base.assoc | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' | awk '{print $",which(head.disc=="SNP"),
					",$",which(head.disc=="A1"),",$",which(head.disc=="BETA"),"}' > rawfile.raw", sep = ""))
				system(paste("awk '{print $",which(head.disc=="SNP"),",$",which(head.disc=="P"),"}' cleaned_base.assoc | tail -n +2 | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' > rangelist_ranges"))
			}	
	        if(length(dir()[grep("cleaned_base.assoc", dir())]) == 0 ){
	            cat("ERROR: Clumping Failed, check dicovery data set format \n "); quit()
	        }
	
	    }
	    if(clump.snps & !dosage & geno.as.list){
			## clump and extract independant snps
			for(chrnum in 1:22){
			    system(paste(plink," --noweb --bfile non_synonymous_snps_only_",chrnum," --clump cleaned_base.assoc --clump-p1 ",clump.p1," --clump-p2 ",clump.p2," --clump-r2 ",clump.r2," --clump-kb ",clump.kb," --out cleaned_base_",chrnum, sep=""))
			}
			system("cat cleaned_base_*clumped > cleaned_base.clumped")
	        system("tail -n +2 cleaned_base.clumped | awk '{print $3}'  | sort -k1,1 | awk '($1 != \"\"){print}'  > LE_SNPs")
			## linkage equilibrium snps
			if(!use.beta){
			    system(paste("tail -n +2 cleaned_base.assoc | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' | awk '{print $",which(head.disc=="SNP"),
				",$",which(head.disc=="A1"),",log($",which(head.disc=="OR"),")}' > rawfile.raw", sep = ""))
			    system(paste("awk '{print $",which(head.disc=="SNP"),",$",which(head.disc=="P"),"}' cleaned_base.assoc | tail -n +2 | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' > rangelist_ranges"))
			}
			if(use.beta){
			    system(paste("tail -n +2 cleaned_base.assoc | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' | awk '{print $",which(head.disc=="SNP"),
				",$",which(head.disc=="A1"),",$",which(head.disc=="BETA"),"}' > rawfile.raw", sep = ""))
			    system(paste("awk '{print $",which(head.disc=="SNP"),",$",which(head.disc=="P"),"}' cleaned_base.assoc | tail -n +2 | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' > rangelist_ranges"))
			}
	        if(length(dir()[grep("cleaned_base.assoc", dir())]) == 0 ){
	            cat("ERROR: Clumping Failed, check dicovery data set format \n "); quit()
	        }
				
	    }
	

		if(clump.snps & dosage & dosage.format == "gen"){
	        ## clump and extract independant snps
	        system(paste(plink," --noweb --dosage ", gen.name, " noheader format=",dos.format," dose",dos.coding," skip0=",dos.skip0," skip1=",dos.skip1," --fam ", target.fam,"  --clump cleaned_base.assoc --clump-p1 ",clump.p1," --clump-p2 ",clump.p2," --clump-r2 ",clump.r2," --clump-kb ",clump.kb," --out cleaned_base", sep=""))
	        system("tail -n +2 cleaned_base.clumped | awk '{print $3}'  | sort -k1,1 | awk '($1 != \"\"){print}'  > LE_SNPs")
	        ## linkage equilibrium snps
	        if(!use.beta){
	            system(paste("tail -n +2 cleaned_base.assoc | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' | awk '{print $",which(head.disc=="SNP"),
	     ",$",which(head.disc=="A1"),",log($",which(head.disc=="OR"),")}'  | sort -k1,1 > rawfile.raw", sep = ""))
	            system(paste("awk '{print $",which(head.disc=="SNP"),",$",which(head.disc=="P"),"}' cleaned_base.assoc | tail -n +2 | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-''  > rangelist_ranges"))
	        }
	        if(use.beta){
	            system(paste("tail -n +2 cleaned_base.assoc | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' | awk '{print $",which(head.disc=="SNP"),
	      ",$",which(head.disc=="A1"),",$",which(head.disc=="BETA"),"}  | sort -k1,1' > rawfile.raw", sep = ""))
	            system(paste("awk '{print $",which(head.disc=="SNP"),",$",which(head.disc=="P"),"}' cleaned_base.assoc | tail -n +2 | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' > rangelist_ranges"))
	        }	
	    }
	
	
	    if(prune.snps){
	        cat(" ################################# \n # \n #   Prune \n # \n ################################# \n")
	        if(!dosage & !(geno.as.list)){
	            ## prune and extract independant snps
	            system(paste(plink," --noweb --bfile non_synonymous_snps_only --indep-pairwise ",prune.kb.wind," kb ",prune.kb.step, prune.kb.r2," --out cleaned_base", sep=" "))
	            system("tail -n +2 cleaned_base.prune.in | awk '{print $1}'  | sort -k1,1 | awk '($1 != \"\"){print $0}'  > LE_SNPs")
	            ## linkage equilibrium snps
	        }
	
	        if(!dosage & (geno.as.list)){
	            ## prune and extract independant snps
	            for(chrnum in 1:22){
	                system(paste(plink," --noweb --bfile non_synonymous_snps_only_",chrnum," --indep-pairwise ",prune.kb.wind," kb ",prune.kb.step," ", prune.kb.r2," --out cleaned_base_",chrnum, sep=""))
	             }
	            system("cat cleaned_base_*prune.in > cleaned_base.prune.in")
	        }
	  
	        system("tail -n +2 cleaned_base.prune.in | awk '{print $1}'  | sort -k1,1 | awk '($1 != \"\"){print $0}'  > LE_SNPs")
	        ## linkage equilibrium snps
	        if(!use.beta){
	            system(paste("tail -n +2 cleaned_base.assoc | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' | awk '{print $",which(head.disc=="SNP"),
		        ",$",which(head.disc=="A1"),",log($",which(head.disc=="OR"),")}' > rawfile.raw", sep = ""))
				system(paste("awk '{print $",which(head.disc=="SNP"),",$",which(head.disc=="P"),"}' cleaned_base.assoc | tail -n +2 | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' > rangelist_ranges"))
            }
	        if(use.beta){
	            system(paste("tail -n +2 cleaned_base.assoc | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' | awk '{print $",which(head.disc=="SNP"),
			    ",$",which(head.disc=="A1"),",$",which(head.disc=="BETA"),"}' > rawfile.raw", sep = ""))
				system(paste("awk '{print $",which(head.disc=="SNP"),",$",which(head.disc=="P"),"}' cleaned_base.assoc | tail -n +2 | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' > rangelist_ranges"))
	        }	
	    }
	
	
	    if(prune.snps & dosage & dosage.format == "gen"){
	        ## prune and extract independant snps
	        system(paste(plink," --noweb --dosage ", gen.name, " noheader format=",dos.format," dose",dos.coding," skip0=",dos.skip0," skip1=",dos.skip1," --fam ", target.fam,"  --indep-pairwise ",prune.kb.wind," ",prune.kb.step," ",prune.kb.r2," --out cleaned_base", sep=""))
	        system("tail -n +2 cleaned_base.clumped | awk '{print $3}'  | sort -k1,1 | awk '($1 != \"\"){print}'  > LE_SNPs")
	        ## linkage equilibrium snps
	        if(!use.beta){
	            system(paste("tail -n +2 cleaned_base.assoc | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' | awk '{print $",which(head.disc=="SNP"),
	            ",$",which(head.disc=="A1"),",log($",which(head.disc=="OR"),")}'  | sort -k1,1 > rawfile.raw", sep = ""))
	            system(paste("awk '{print $",which(head.disc=="SNP"),",$",which(head.disc=="P"),"}' cleaned_base.assoc | tail -n +2 | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-''  > rangelist_ranges"))
	        }
	        if(use.beta){
	            system(paste("tail -n +2 cleaned_base.assoc | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' | awk '{print $",which(head.disc=="SNP"),
	            ",$",which(head.disc=="A1"),",$",which(head.disc=="BETA"),"}  | sort -k1,1' > rawfile.raw", sep = ""))
	            system(paste("awk '{print $",which(head.disc=="SNP"),",$",which(head.disc=="P"),"}' cleaned_base.assoc | tail -n +2 | sort -k1,1 | join -1 1 -2 1 LE_SNPs ``-'' > rangelist_ranges"))
	        }	
	    }
	
	
	
	    if(!clump.snps & !prune.snps){
		    if(!use.beta){
			    system(paste("tail -n +2 cleaned_base.assoc | awk '{print $",which(head.disc=="SNP"),
				",$",which(head.disc=="A1"),",log($",which(head.disc=="OR"),")}' | sort -k1,1 > rawfile.raw", sep = ""))
			    system(paste("awk '{print $",which(head.disc=="SNP"),",$",which(head.disc=="P"),"}' cleaned_base.assoc | tail -n +2   > rangelist_ranges"))
		    }
		    if(use.beta){
			    system(paste("tail -n +2 cleaned_base.assoc | awk '{print $",which(head.disc=="SNP"),
				",$",which(head.disc=="A1"),",$",which(head.disc=="BETA"),"}' | sort -k1,1 > rawfile.raw", sep = ""))
			    system(paste("awk '{print $",which(head.disc=="SNP"),",$",which(head.disc=="P"),"}' cleaned_base.assoc | tail -n +2  > rangelist_ranges"))
		    }
	    }
	
	
	
	    cat(" ################################# \n # \n #   Deal with strand flips if target is in genotype format and produce input files for polygenic scoring \n # \n ################################# \n")
	
	    if(!dosage){
	        ## strand flips
	        system(paste("awk '{print $2,$5,$6}' non_synonymous_snps_only.bim | sort -k1,1 | join -1 1 -2 1   rawfile.raw   ``-'' | awk '($2 != $4 && $2 != $5) {print $1}' > flip_list.txt", sep = 	""))
            ## Complete records
	        system(paste("awk '{print $2}' non_synonymous_snps_only.bim | sort -k1,1 | join -1 1 -2 1   rangelist_ranges   ``-''  > Complete_Allele_List.txt", sep = 	""))
	        ComAlLi <- read.table("Complete_Allele_List.txt", head=F)
            names(ComAlLi) <- c("SNP", "P")
	  	  	  
	        if(!(geno.as.list)){ 
	            system(paste( plink, "        --noweb --bfile non_synonymous_snps_only --flip flip_list.txt ",mhc,"  --make-bed   --out flipped_target", sep = ""))
	        }
	        if(geno.as.list){
	  	        for(chrnum in 1:5){ 
	                system(paste( plink, "        --noweb --bfile non_synonymous_snps_only_",chrnum," --flip flip_list.txt  --make-bed   --out flipped_target_",chrnum, sep = ""))
	            }
	  	        for(chrnum in 6){ 
	                system(paste( plink, "        --noweb --bfile non_synonymous_snps_only_",chrnum," --flip flip_list.txt ",mhc,"  --make-bed   --out flipped_target_",chrnum, sep = ""))
	            }
	  	        for(chrnum in 7:22){ 
	                system(paste( plink, "        --noweb --bfile non_synonymous_snps_only_",chrnum," --flip flip_list.txt   --make-bed   --out flipped_target_",chrnum, sep = ""))
	            }
	        }
	
	        system("rm rangelist.txt", ignore.stdout=T,ignore.stderr=T)
	        if(!fastscore){
	            system( capture.output(cat("awk 'BEGIN{ for (i=",slower,"; i < ",supper,"; i+=",sinc,") printf(\"%.",ceiling(-log10(sinc)),"f\\n\", i); }' | awk '{print $1,0,$1}' > 	rangelist.txt \n ", sep=""))[1])
	            lists <- read.table("rangelist.txt", head = F)
	            system( capture.output(cat("awk 'BEGIN{ for (i=",slower,"; i < ",supper,"; i+=",sinc,") printf(\"%.",ceiling(-log10(sinc)),"f\\n\", i); }' |  awk '{OFS=\".\"}			{print \"PROFILES\",$1,\"profile\"}' > profile_list 	\n ", sep = ""))[1])
	        }
	        if(fastscore){
	            write.table(data.frame(barchart.levels, rep(0, times=length(barchart.levels)), barchart.levels), "rangelist.txt", col.names = F, row.names = F, quote = F)
	            lists <- read.table("rangelist.txt", head = F)
	            write.table(data.frame(paste("PROFILES",barchart.levels,"profile",sep=".")), "profile_list", col.names = F, row.names = F, quote = F)
	        } 
	        if(mend.score){
	            system(paste("tail -n +2 cleaned_base.assoc | sort -gk",which(head.disc=="P"),",",which(head.disc=="P")," | head -n ", mend.score.len," | awk      '{print $",which(head.disc=="P"),"}' > best.n.snps ", sep=""))
	            best.n <- read.table("best.n.snps",head=F)
	            best.n <- data.frame(best.n$V1, 0, best.n$V1)
	            write.table(best.n, "best.n", col.names=F, row.names=F, quote=F)
	            system("mv rangelist.txt temp.rangelist.txt")
	            system("cat temp.rangelist.txt best.n > rangelist.txt")
	            lists <- read.table("rangelist.txt", head = F)
	            system("mv profile_list temp.1.profile_list")
	            system(paste("awk '{print $1}' best.n.snps | awk '{OFS=\".\"}	{print \"PROFILES\",$1,\"profile\"}' >  temp.2.profile_list"))
	            system("cat temp.1.profile_list temp.2.profile_list > profile_list" )
	        }
	        if(score.at.1){
	            system("mv profile_list temp.profile_list")
	            system("cat temp.profile_list PROFILES.1.profile > profile_list")
	            system("mv rangelist.txt temp.rangelist.txt")
	            system("echo '1 0 1' | cat temp.rangelist.txt ``-'' > rangelist.txt")
            } 
        }  
	  
	    if(dosage){
	        system("rm rangelist.txt", ignore.stdout=T,ignore.stderr=T)
	        if(!fastscore){
	            system( capture.output(cat("awk 'BEGIN{ for (i=",slower,"; i < ",supper,"; i+=",sinc,") printf(\"%.",ceiling(-log10(sinc)),"f\\n\", i); }' | awk '{print $1,0,$1}' > 	  rangelist.txt \n ", sep=""))[1])
	            lists <- read.table("rangelist.txt", head = F)
	            write.table(data.frame(paste("PROFILES.S", seq(1, length(lists$V1), 1), ".profile", sep = "")), "profile_list", col.names=F, row.names=F,quote=F)	
	        }
	        if(fastscore){
	            write.table(data.frame(barchart.levels, rep(0, times=length(barchart.levels)), barchart.levels), "rangelist.txt", col.names = F, row.names = F, quote = F)
	            lists <- read.table("rangelist.txt", head = F)
	            write.table(data.frame(paste("PROFILES.S",seq(1, length(lists$V1), 1),".profile",sep="")), "profile_list", col.names = F, row.names = F, quote = F)
	        }
	        if(mend.score){
	            system(paste("tail -n +2 cleaned_base.assoc | sort -gk",which(head.disc=="P"),",",which(head.disc=="P")," | head -n ", mend.score.len," | awk      '{print $",which(head.disc=="P"),"}' > best.n.snps ", sep=""))
	            best.n <- read.table("best.n.snps",head=F)
	            best.n <- data.frame(best.n$V1, 0, best.n$V1)
	            write.table(best.n, "best.n", col.names=F, row.names=F, quote=F)
	            system("mv rangelist.txt temp.rangelist.txt")
	            system("cat temp.rangelist.txt best.n > rangelist.txt")
	            lists <- read.table("rangelist.txt", head = F)
	            system("mv profile_list temp.1.profile_list")
	            write.table(data.frame(paste("PROFILES.S", seq(length(lists$V1)+1, length(lists$V1)+mend.score.len, 1), ".profile", sep = "")), "temp.2.profile_list", col.names=F, row.names=F,quote=F)	
	            system("cat temp.1.profile_list temp.2.profile_list > profile_list" )
	        }
	        if(score.at.1){
	            system("mv profile_list temp.profile_list")
        	    if(!mend.score){
	                system("cat temp.profile_list PROFILES.S",length(lists$V1)+mend.score.len+1,".profile > profile_list")
	            }
	            if(!mend.score){
	                system("cat temp.profile_list PROFILES.S",length(lists$V1)+1,".profile > profile_list")
	            }   
	            system("mv rangelist.txt temp.rangelist.txt")
	            system("echo '1 0 1' | cat temp.rangelist.txt ``-'' > rangelist.txt")
	        }
	    }
	
	    cat(" ################################# \n # \n #   Polygenic scoring! \n # \n ################################# \n")
	
	
	    if(!dosage & !(geno.as.list)){
		    system(paste(plink, "        --noweb --bfile flipped_target   --score rawfile.raw   --q-score-range rangelist.txt rangelist_ranges  --out PROFILES"))
	    }
	
	    if(!dosage & (geno.as.list)){
	        for(chrnum in 1:22){
	            system(paste(plink, "        --noweb --bfile flipped_target_",chrnum,"   --score rawfile.raw   --q-score-range rangelist.txt rangelist_ranges  --out ",chrnum,"PROFILES",sep=""))
	        }
	    }
        if(dosage & is.na(dos.path.to.lists) & is.na(dos.list.file) & dos.coding == 1){
            system(paste(plink," --noweb --dosage ", gen.name, " noheader format=",dos.format," dose",dos.coding," skip0=",dos.skip0," skip1=",dos.skip1," --fam ", target.fam,"  --score rawfile.raw --q-score-file rangelist_ranges --q-score-range rangelist.txt   --out PROFILES", sep = ""))
        }
        if(dosage & is.na(dos.path.to.lists) & is.na(dos.list.file) & dos.coding == 2){
            system(paste(plink," --noweb --dosage ", gen.name, " noheader format=",dos.format," skip0=",dos.skip0," skip1=",dos.skip1," --fam ", target.fam,"  --score rawfile.raw --q-score-file rangelist_ranges --q-score-range rangelist.txt   --out PROFILES", sep = ""))
        }
        if(dosage & is.na(dos.path.to.lists) &is.na(dos.list.file) & dos.coding != 1 & dos.coding != 2 ){
            system(paste(plink," --noweb --dosage ",gen.name, " noheader format=",dos.format," dose",dos.coding," skip0=",dos.skip0," skip1=",dos.skip1," --fam ", target.fam,"  --score rawfile.raw --q-score-file rangelist_ranges --q-score-range rangelist.txt   --out PROFILES", sep = ""))
        }

	    ## NB: score chromosome by chromosome if dosage data in separate files per chromosome. Sum these later on to get profiles
	    if(dosage & (!is.na(dos.path.to.lists) | !is.na(dos.list.file))){
	        for(i in 1:22){
	            if(!is.na(dos.path.to.lists)){
	                gen.name <- paste(gsub(CHRNAME, i, dos.path.to.lists))
	            }
	            if(!is.na(dos.list.file)){
	                gen.name <- as.character(read.table(dos.list.file, head = F)[i,1])
	            }
	            system(paste(plink," --noweb --dosage ",gen.name," noheader format=",dos.format," dose",dos.coding," skip0=",dos.skip0," skip1=",dos.skip1," --fam ", target.fam,"  --score rawfile.raw --q-score-file rangelist_ranges --q-score-range rangelist.txt   --out ",i,"PROFILES", sep = ""))  
	        }
	    }            
	
	
	    if(covary & is.na(user.covariate.file)){
	        cat(" ################################# \n # \n #   Covary by generated dimensions \n # \n ################################# \n")
	    }
	
	    if(covary & is.na(user.covariate.file)){
	        if(!dosage){
	  	        if(geno.as.list){
	  		        cat("ERROR: Please merge genome to calculate ancenstry-informative dimensions \n");quit()
	  	        }
	  	        if(!(geno.as.list)){
                    system("echo 6 26000000 33000000 mhc > mhc.txt")
                    system(paste(plink, " --noweb --bfile ",target, "         --exclude range  mhc.txt   --make-bed  --out target_no_mhc "))
                    system(paste(plink, " --noweb --bfile   target_no_mhc     --indep-pairwise 100 25 0.2     --out prune_target"))
                    system(paste(plink, " --noweb --bfile ",target, "         --extract prune_target.prune.in    --make-bed   --out prune_target_out"))
                    system(paste(plink, " --noweb --bfile   prune_target_out  --genome     --out prune_target_out"))
                    if(ancestry.dim == "MDS"){
	                    system(paste(plink, " --noweb --bfile   prune_target_out  --read-genome prune_target_out.genome   --cluster --mds-plot", length(covariates), " --out 	  ANCESTRY_INFORMATIVE_DIMENSIONS"))
	                    pc.file <- "ANCESTRY_INFORMATIVE_DIMENSIONS.mds"
	                }
	                if(ancestry.dim == "PCA"){
	                    system(paste(plink, " --noweb --bfile   prune_target_out  --read-genome prune_target_out.genome  --cluster --pca ", length(covariates), " header --out 		ANCESTRY_INFORMATIVE_DIMENSIONS"))
	                    system("sed -e '1s/P//g' ANCESTRY_INFORMATIVE_DIMENSIONS.eigenvec > ANCESTRY_INFORMATIVE_DIMENSIONS.readable")
	                    pc.file <- "ANCESTRY_INFORMATIVE_DIMENSIONS.readable"
	                }
	            }
	        }
	    }
	
	    if(covary & !is.na(user.covariate.file)){
	        pc.file<- user.covariate.file
	    }
	
	    if(ext.phen & !multiple.target.phenotypes){
	        pheno.data <- read.table(pheno.file, head = F)
	    }
	    if(ext.phen & multiple.target.phenotypes){
	        pheno.data <- read.table(pheno.file, head = T)
	    }
	
	    ####################################################
	    ## Extract phenotype data and save it somewhere::
        ####################################################
        if(quantiles & dosage & is.na(ext.phen)){
            cat("WARNING: Cannot produce quantiles plot for dosage data \n unless external phenotype file is supplied \n")
            quantiles <- F
        }
	    if(quantiles){
	        if(ext.phen){
	            phen.file.internal <- pheno.data 
	        }
	        if(!ext.phen){
	            if(!dosage){
	                if(!geno.as.list){
                        phen.file.internal <- read.table(paste(target, "fam", sep  = "."), head = F)[,c("V2", "V6")] 
	                }
	                if(geno.as.list){
                        phen.file.internal <- read.table("flipped_target_1.fam", head = F)[,c("V2", "V6")]     	
	                }
	            }
            }
            names(phen.file.internal) <- c("ID", "PHEN")
            if(levels(as.factor(phen.file.internal$PHEN))[1]=="1"&levels(as.factor(phen.file.internal$PHEN))[2]=="2"){ 
                phen.file.internal$PHEN <- phen.file.internal$PHEN - 1
	        }  
            phen.file.internal <- phen.file.internal[phen.file.internal$PHEN != -9 , ]
        }
	
	    if(covary){	
	        pcs <- read.table(pc.file, head = T)
	    }
	    profile.list <- read.table("profile_list", head = F)
	    lists.full <- lists	
		
	    ## SOMETIMES: an empty list will not print. this section is important - it takes the lists which did print and retains only these from the list of names to read and test model fit etc
	    if(!dosage){
	        if(geno.as.list){
	            for(i in 1:22){
	                if(i == 1){
	                    lists.temp <- lists[paste(i,profile.list$V1,sep="") %in% list.files(pattern='*PROFILES.*\\.profile'),]
	                    lists.out <- lists.temp
	                }
	                if(i > 1){
	                    lists.temp <- lists[paste(i,profile.list$V1,sep="") %in% list.files(pattern='*PROFILES.*\\.profile'),]
	                    lists.out <- rbind(lists.out, lists.temp)
	                }
                }
	            lists <- lists[!duplicated(lists$V1),]  
	            profile.list.temp <- data.frame(profile.list, lists.full)
	            names(profile.list.temp) <- c("V1", "V2", "V3","V4")
	            reduced.list <- data.frame(profile.list.temp$V1[profile.list.temp$V2 %in% lists$V1])
	            names(reduced.list) <- c("V1")
	        }
	        if(!geno.as.list){
        	    lists <- lists[profile.list$V1 %in% list.files(pattern='PROFILES.*\\.profile'),]
	            reduced.list <- profile.list$V1[profile.list$V1 %in% list.files(pattern='PROFILES.*\\.profile')]
	        }
        }
	    if(dosage){
		    lists <- lists[profile.list$V1 %in% list.files(pattern='PROFILES.*\\.profile'),]
		    reduced.list <- profile.list$V1[profile.list$V1 %in% list.files(pattern='PROFILES.*\\.profile')]
	    }
	    if(!no.regression){
	        cat(" ################################# \n # \n #   Regression Models \n # \n ################################# \n")
	    }

	    ranges.list <- read.table("rangelist_ranges",head=F)

	    if(!multiple.target.phenotypes){
	        p.out <- as.vector(1)
	        r2.out <- as.vector(1)
	        nsnps <- as.vector(1)
	        coefficient <- as.vector(1)
	        s.err <- as.vector(1)
	        prof.temp <- as.vector(1)
	        prev.files <- F
	        abc.r2 <- as.vector(1)
	        ## logistic regression on phenotype
	        if(binary.target){
        	    for(i in 1:length(lists$V1)){
	                if(dosage & (!is.na(dos.path.to.lists) | !is.na(dos.list.file))){
	                    if(!is.na(dos.list.file)){
	                        max.files <- dim(read.table(dos.list.file))[1]
	                    }
	                    if(!is.na(dos.path.to.lists)){
	                        max.files <- 22
	                    }
	                    for(j in 1:max.files){
	                        if(paste(j, profile.list[i,1], sep = "") %in% list.files()){
	                            prof.temp <- read.table(paste(j, profile.list[i,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	                            names(prof.temp) <- c("temp.IID", "temp.SCORE","temp.PHENO")
	                            if(prev.files){
	                                prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	                                prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	                                prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	                            }
	                            prof.out <- prof.temp
	                            names(prof.out) <- c("IID","SCORE","PHENO")
	                            prof <- prof.out
	                            prev.files <- T
	                        }
	                    }
	                }
	                if(geno.as.list){
	                    for(j in 1:22){
	                        if(paste(j, profile.list[i,1], sep = "") %in% list.files()){
	                            prof.temp <- read.table(paste(j, profile.list[i,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	                            names(prof.temp) <- c("temp.IID", "temp.SCORE","temp.PHENO")
	                            if(prev.files){
	                                prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	                                prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	                                prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	                            }
	                            prof.out <- prof.temp
	                            names(prof.out) <- c("IID","SCORE","PHENO")
	                            prof <- prof.out
	                            prev.files <- T         
	                        }
	                        prev.files <- F
	                    }
	                }
	                if((!dosage & !geno.as.list)  | (dosage & is.na(dos.path.to.lists) & is.na(dos.list.file))){			
		                ## Default - no fastscore
		                prof <- read.table(paste(reduced.list[i], sep = ""), head  =T)
	                    prof <- prof[,c("IID", "PHENO", "SCORE")]
	                }
	                if(ext.phen){
	                    prof <- subset(prof, select=-c(PHENO))
		                prof <- merge(x = prof, by.x = "IID", y = pheno.data, by.y = "V1")
	                }
		            if(!ext.phen) {
	                    prof <- prof[prof$PHENO != -9,] 
	                    prof <- prof[!is.na(prof$PHENO),]
	                    if(levels(as.factor(prof$PHENO))[1]=="1"&levels(as.factor(prof$PHENO))[2]=="2"){ 
		                    prof$PHENO <- prof$PHENO - 1
	                    }
	                    if(!(levels(as.factor(prof$PHENO))[1]=="0"&levels(as.factor(prof$PHENO))[2]=="1")){ 
		                    print("ERROR: Unrecognised values for binary phenotype.");quit()
	                    }
	                }
	                if(length(levels(as.factor(prof$PHENO)) )< 2 & !ext.phen){
	                    cat("ERROR: Phenotype does not have more than one level. Will not perform regression \n")
	                    no.regression <- T
	                }
	                if(!no.regression){
	                    if(covary & ext.phen){
	                        prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
	                        model.logit <- glm(V2 ~., family="binomial", data = prof[,c("V2","SCORE", covariates)])
	                        model.null <-  glm(V2 ~., family="binomial", data = prof[,c("V2", covariates)])
	                        p.out[i]  <- summary(model.logit)$coefficients[2,4]
	                        if(for.meta){
	                            coefficient[i]  <- summary(model.logit)$coefficients[2,1]
	                            s.err[i] <- summary(model.logit)$coefficients[2,2]
	                        }
		                    r2.out[i] <- NagelkerkeR2(model.logit)$R2 - NagelkerkeR2(model.null)$R2  
	                    }
	                    if(covary & !ext.phen){
	                        prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
		                    model.logit <- glm(PHENO ~ ., family="binomial", data = prof[,c("PHENO", "SCORE", covariates)])
		                    model.null <-  glm(PHENO  ~ ., family="binomial", data = prof[,c("PHENO", covariates)])
		                    p.out[i]  <- summary(model.logit)$coefficients[2,4]
	                        if(for.meta){
	                            coefficient[i]  <- summary(model.logit)$coefficients[2,1]
	                            s.err[i] <- summary(model.logit)$coefficients[2,2]
	                        }
		                    r2.out[i] <- NagelkerkeR2(model.logit)$R2 - NagelkerkeR2(model.null)$R2
	                    }
	                    if(!covary & ext.phen){
	                        p.out[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="binomial")))$coefficients[2,4]
	                        if(for.meta){
	                            coefficient[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="binomial")))$coefficients[2,1]
	                            s.err[i] <- summary(with(prof, glm(V2  ~ SCORE, family="binomial")))$coefficients[2,2]
	                        }
	                        r2.out[i] <- NagelkerkeR2(with(prof, glm(V2 ~ SCORE, family="binomial")))$R2
	                    }
	                    if(!covary & !ext.phen){
	                        p.out[i]  <- summary(with(prof, glm(PHENO  ~ SCORE, family="binomial")))$coefficients[2,4]
	                        if(for.meta){
	                            coefficient[i]  <- summary(with(prof, glm(PHENO  ~ SCORE, family="binomial")))$coefficients[2,1]
	                            s.err[i] <- summary(with(prof, glm(PHENO  ~ SCORE, family="binomial")))$coefficients[2,2]
	                        }
	                        r2.out[i] <- NagelkerkeR2(with(prof, glm(PHENO ~ SCORE, family="binomial")))$R2
	                    }
	                }
	                prof.red <- prof[,c("IID","SCORE")]
	                names(prof.red) <- c("IID", paste("pT", lists$V1[i], sep = "_"))
	                if(report.individual.scores){
	                    if(i == 1){
	                        prof.all.scores <- prof.red
	                    }
	                    if(i > 1){
	                        prof.all.scores <- merge(x = prof.all.scores, by.x = "IID", y = prof.red, by.y = "IID")    
	                    }
	                }
	                nsnps[i] <- dim(ranges.list[ranges.list$V2 < lists$V1[i],])[1]
	                if(i %in% (seq(1, length(lists$V1), 1)[seq(1, length(lists$V1), 1) %in% ceiling(length(lists$V1)*seq(0, 1,0.1))] )){
	                    cat(paste("Regression Models: " ,round(100*(i/length(lists$V1)),digits=-1), "% Complete \n",sep=""))
	                }
                    if(calculate.abc.r2){
          	            abc.r2[i] <- estimateVg2FromP(prevalence1=prevalence1,
          	                                            prevalence2=prevalence2,
          	                                            n1=n1, 
          	                                            sampling1= sampling1, 
          	                                            n2= n2, 
          	                                            sampling2= sampling2, 
          	                                            nsnp= nsnps[i], 
          	                                            plower = 0, 
          	                                            pupper = lists$V1[i], 
          	                                            binary=binary.target, 
          	                                            p=p.out[i],
          	                                            nullfraction=pi0
          	                                            )$vg
                    }
	            }
	        }	
	        ## linear regression on phenotype
	        if(!binary.target){
	            for(i in 1:length(lists$V1)){
	                if(dosage & (!is.na(dos.path.to.lists) | !is.na(dos.list.file))){
	                    for(j in 1:22){
	                        prof.temp <- read.table(paste(j, profile.list[j,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	                        names(prof.temp) <- c("temp.IID", "temp.SCORE")
	                        if(j > 1){
	                            prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	                            prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	                            prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	                        }
	                        prof.out <- prof.temp
	                        names(prof.out) <- c("IID","SCORE","PHENO")
	                    }
    	                prof <- prof.out
	                }
	                if(!dosage | (is.na(dos.path.to.lists) & is.na(dos.list.file))){			
	                    ## Default - no fastscore
	                    prof <- read.table(paste(reduced.list[i], sep = ""), head  =T)
	                    prof <- prof[,c("IID", "PHENO", "SCORE")]
	                }
	                if(ext.phen){
	                    prof <- subset(prof, select=-c(PHENO))
	                    prof <- merge(x = prof, by.x = "IID", y = pheno.data, by.y = "V1")
	                } 
	                if(length(levels(as.factor(prof$PHENO)) )< 2 & !ext.phen){
	                    cat("ERROR: Phenotype does not have more than one level. Will not perform regression \n")
	                    no.regression <- T
	                }
	                if(!no.regression){
	                    if(covary & ext.phen){
	                        prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
		                    model.logit <- glm(V2  ~ ., family="gaussian", data = prof[,c("V2", "SCORE", covariates)])
		                    model.null <-  glm(V2  ~ ., family="gaussian", data = prof[,c("V2", covariates)])
		                    p.out[i]  <- summary(model.logit)$coefficients[2,4]
	                        if(for.meta){
	                            coefficient[i]  <- summary(model.logit)$coefficients[2,1]
	                            s.err[i] <- summary(model.logit)$coefficients[2,2]
	                        }
		                    r2.out[i] <- (1 - (sum( (prof$V2[complete.cases(prof[,c("V2", "SCORE", covariates)])]-predict(model.logit))^2 ) / sum( (prof$V2[complete.cases(prof[,c("V2", "SCORE", covariates)])]-mean(prof$V2[complete.cases(prof[,c("V2", "SCORE", covariates)])]))^2 ) ) )  - (1 - ( sum( (prof$V2[complete.cases(prof[,c("V2", "SCORE", covariates)])]-predict(model.null))^2 ) / sum( (prof$V2[complete.cases(prof[,c("V2", "SCORE", covariates)])]-mean(prof$V2[complete.cases(prof[,c("V2", "SCORE", covariates)])]))^2 ) ) )  
		                }
		                if(covary & !ext.phen){
	                        prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
	                        model.logit <- glm(PHENO ~ ., family="gaussian", data = prof[,c("PHENO", "SCORE", covariates)])
	                        model.null <-  glm(PHENO ~ ., family="gaussian", data = prof[,c("PHENO", covariates)])
	                        p.out[i]  <- summary(model.logit)$coefficients[2,4]
	                        if(for.meta){
	                            coefficient[i]  <- summary(model.logit)$coefficients[2,1]
	                            s.err[i] <- summary(model.logit)$coefficients[2,2]
	                        }
	                        r2.out[i] <- (1 - ( sum( (prof$PHENO[complete.cases(prof[,c("PHENO", "SCORE", covariates)])] - 
	                                                predict(model.logit))^2 ) / 
	                                        sum( (prof$PHENO[complete.cases(prof[,c("PHENO", "SCORE", covariates)])] -
	                                            mean(prof$PHENO[complete.cases(prof[,c("PHENO", "SCORE", covariates)])]))^2 ) ) )  - 
	                                    (1 - ( sum( (prof$PHENO[complete.cases(prof[,c("PHENO", "SCORE", covariates)])] - 
	                                                 predict(model.null))^2 ) / 
	                                        sum( (prof$PHENO[complete.cases(prof[,c("PHENO", "SCORE", covariates)])] - 
	                                            mean(prof$PHENO[complete.cases(prof[,c("PHENO", "SCORE", covariates)])]))^2 ) ) ) 
	                    }
	                    if(!covary & ext.phen){
	                        p.out[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="gaussian")))$coefficients[2,4]
	                        if(for.meta){
	                            coefficient[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="gaussian")))$coefficients[2,1]
	                            s.err[i] <- summary(with(prof, glm(V2  ~ SCORE, family="gaussian")))$coefficients[2,2]
	                        }
	                        r2.out[i] <-  1 - (sum( (prof$V2[complete.cases(prof[,c("V2", "SCORE")])]-predict(with(prof[complete.cases(prof[,c("V2", "SCORE")]),], glm(V2 ~ SCORE, family="gaussian"))))^2 )/sum( (prof$V2[complete.cases(prof[,c("V2", "SCORE")])]-mean(prof$V2[complete.cases(prof[,c("V2", "SCORE")])]))^2))   			
	                    }
	                    if(!covary & !ext.phen){
	                        p.out[i]  <- summary(with(prof, glm(PHENO  ~ SCORE, family="gaussian")))$coefficients[2,4]
	                        if(for.meta){
	                            coefficient[i]  <- summary(with(prof, glm(PHENO  ~ SCORE, family="gaussian")))$coefficients[2,1]
	                            s.err[i] <- summary(with(prof, glm(PHENO  ~ SCORE, family="gaussian")))$coefficients[2,2]
	                        }
	                        r2.out[i] <- 1 - (sum( (prof$PHENO[complete.cases(prof[,c("PHENO", "SCORE")])]-predict(with(prof[complete.cases(prof[,c("PHENO", "SCORE")]),], glm(PHENO ~ SCORE, family="gaussian"))))^2 )/sum( (prof$PHENO[complete.cases(prof[,c("PHENO", "SCORE")])]-mean(prof$PHENO[complete.cases(prof[,c("PHENO", "SCORE")])]))^2 )) 
	                    }
	                }
	                prof.red <- prof[,c("IID","SCORE")]
	                names(prof.red) <- c("IID", paste("pT", lists$V1[i], sep = "_"))
	                if(report.individual.scores){
	                    if(i == 1){
	                        prof.all.scores <- prof.red
	                    }
	                    if(i > 1){
	                        prof.all.scores <- merge(x = prof.all.scores, by.x = "IID", y = prof.red, by.y = "IID")    
	                    }
	                }
	                nsnps[i] <- dim(ranges.list[ranges.list$V2 < lists$V1[i],])[1]
	                if(i %in% (seq(1, length(lists$V1), 1)[seq(1, length(lists$V1), 1) %in% ceiling(length(lists$V1)*seq(0, 1,0.1))] )){
	                    cat(paste("Regression Models: " ,round(100*(i/length(lists$V1)),digits=-1), "% Complete \n",sep=""))
	                }
                    if(calculate.abc.r2){
          	            abc.r2[i] <- estimateVg2FromP(prevalence1=prevalence1,
                                              	        prevalence2=prevalence2,
                                              	        n1=n1, 
                                              	        sampling1= sampling1, 
                                              	        n2= n2, 
                                              	        sampling2= sampling2, 
                                              	        nsnp= nsnps[i], 
                                              	        plower = 0, 
                                              	        pupper = lists$V1[i], 
                                              	        binary=binary.target, 
                                              	        p=p.out[i],
                                              	        nullfraction=pi0
          	                                        )$vg
                    }
	      
	            }
	        }
	    }

  	    if(emp.pval){
            cat(" ################################# \n # \n #   Calculate Empirical P-value for top threshold \n # \n ################################# \n")
            p.perm <- as.vector(1)
            i <- which.min(p.out)
            if(dosage & (!is.na(dos.path.to.lists) | !is.na(dos.list.file))){
                if(!is.na(dos.list.file)){
                    max.files <- dim(read.table(dos.list.file))[1]
	            }
	            if(!is.na(dos.path.to.lists)){
	                max.files <- 22
	            }
    	        for(j in 1:max.files){
	                if(paste(j, profile.list[i,1], sep = "") %in% list.files()){
	                    prof.temp <- read.table(paste(j, profile.list[i,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	                    names(prof.temp) <- c("temp.IID", "temp.SCORE","temp.PHENO")
	                    if(prev.files){
	                        prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	                        prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	                        prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	                    }
	                    prof.out <- prof.temp
	                    names(prof.out) <- c("IID","SCORE","PHENO")
	                    prof <- prof.out
	                    prev.files <- T
	                }
	            }
	        }
	        if(geno.as.list){
	            for(j in 1:22){
	                if(paste(j, profile.list[i,1], sep = "") %in% list.files()){
	                    prof.temp <- read.table(paste(j, profile.list[i,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	                    names(prof.temp) <- c("temp.IID", "temp.SCORE","temp.PHENO")
	                    if(prev.files){
	                        prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	                        prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	                        prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	                    }
	                    prof.out <- prof.temp
	                    names(prof.out) <- c("IID","SCORE","PHENO")
	                    prof <- prof.out
	                    prev.files <- T         
	                }
	                prev.files <- F
	            }
	        }
	        if((!dosage & !geno.as.list)  | (dosage & is.na(dos.path.to.lists) & is.na(dos.list.file))){			
		        ## Default - no fastscore
		        prof <- read.table(paste(reduced.list[i], sep = ""), head  =T)
	            prof <- prof[,c("IID", "PHENO", "SCORE")]
            }
	        if(ext.phen){
	            prof <- subset(prof, select=-c(PHENO))
		        prof <- merge(x = prof, by.x = "IID", y = pheno.data, by.y = "V1")
	        }
		    if(!ext.phen & binary.target) {
	            prof <- prof[prof$PHENO != -9,] 
	            prof <- prof[!is.na(prof$PHENO),]
	            if(levels(as.factor(prof$PHENO))[1]=="1"&levels(as.factor(prof$PHENO))[2]=="2"){ 
		            prof$PHENO <- prof$PHENO - 1
	            }
	            if(!(levels(as.factor(prof$PHENO))[1]=="0"&levels(as.factor(prof$PHENO))[2]=="1")){ 
		            print("ERROR: Unrecognised values for binary phenotype.");quit()
	            }
	        }

            if(binary.target){
          	    n.binom = dim(prof)[1]
          	    prob.binom = sum(prof$PHENO)/dim(prof)[1]
            }
            if(!binary.target){
          	    mean.norm = mean(prof$PHENO)
            	sd.norm = sd(prof$PHENO)
            }
            for(i in 1:n.emp.perms){
          	    if(binary.target){
                    prof$PHENO <- rbinom(n = n.binom, size = 1, prob = prob.binom)
	                if(covary){
	          	        if(i == 1){
	                        prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
	                    }
	                    model.logit <- glm(PHENO ~., family="binomial", data = prof[,c("PHENO","SCORE", covariates)])
	                    p.perm[i]  <- summary(model.logit)$coefficients[2,4]
	                }
	                if(!covary){
	                    p.perm[i]  <- summary(with(prof, glm(PHENO  ~ SCORE, family="binomial")))$coefficients[2,4]
	                }
          	    }
          	    if(!binary.target){
          		    prof$PHENO <- rnorm(n = n.binom, mean = mean.norm, sd = sd.norm)
	                if(covary){
	                    prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
	                    model.logit <- glm(PHENO ~., family="gaussian", data = prof[,c("PHENO","SCORE", covariates)])
	                    p.perm[i]  <- summary(model.logit)$coefficients[2,4]
	                }
	                if(!covary){
	                    p.perm[i]  <- summary(with(prof, glm(PHENO  ~ SCORE, family="gaussian")))$coefficients[2,4]
	                }
          	    }
                if(report.perm.phen){
          	        if(i == 1){
          	  	        perm.phen.db <- data.frame(prof$IID, prof$PHENO)
          	        }  
          	        if(i > 1){
          	  	        perm.phen.db <- data.frame(perm.phen.db, prof$PHENO)
          	  	        names(perm.phen.db) <- c("IID", paste("PERM", seq(1, i), sep = "_"))
          	        }  
                }
            }
            if(report.perm.phen){
                write.table(perm.phen.db, paste(figname, "PERM_PHENOTYPES.txt", sep="_"), col.names = T, row.names = F, quote = F)
            }
            if(report.perm.pvals){
          	    write.table(data.frame(PERM = paste("PERM", seq(1, i), sep = "_"), PVAL = p.perm), paste(figname, "PERM_PVALS.txt", sep="_"), col.names = T, row.names = F, quote = F)
            }
            if(!emp.alpha){
    	        cat(paste(" # \n # \n # \n # \n # \n #   Empirical P-value at top threshold - pT = ", lists$V1[which.min(p.out)]," - is ", sum(p.perm < p.out[which.min(p.out)]) / length(p.perm) ," \n # \n # \n # \n # \n # \n #\n"))
            }
            if(emp.alpha){
    	        cat(paste(" # \n # \n # \n # \n # \n #   Empirical alpha threshold is ", sort(p.perm)[round(n.emp.perms*0.05)] ," \n # \n # \n # \n # \n # \n #\n"))
            }
  	    }
	    ## Alternately, regression with multiple phenotypes::
	    if(multiple.target.phenotypes){
	        top.thresh <- as.vector(1)
	        top.thresh.pval <- as.vector(1)
            top.thresh.r2 <- as.vector(1)
            top.thresh.r2.dir <- as.vector(1)
            for(k in 1:length(target.phenotypes)){
	  	        cat(paste("Polygenic Scores Regressed on ", target.phenotypes[k] ," Status \n"))
	            binary.target <- target.phenotypes.binary[k]
	            p.out <- as.vector(1)
	            r2.out <- as.vector(1)
	            nsnps <- as.vector(1)
	            coefficient <- as.vector(1)
	            s.err <- as.vector(1)
	            prof.temp <- as.vector(1)
	            prev.files <- F
	            ## logistic regression on phenotype
	            if(binary.target){
	                for(i in 1:length(lists$V1)){
	                    if(dosage & (!is.na(dos.path.to.lists) | !is.na(dos.list.file))){
	                        if(!is.na(dos.list.file)){
	                            max.files <- dim(read.table(dos.list.file))[1]
	                        }
	                        if(!is.na(dos.path.to.lists)){
	                            max.files <- 22
	                        }
	                        for(j in 1:max.files){
	                            if(paste(j, profile.list[i,1], sep = "") %in% list.files()){
	                                prof.temp <- read.table(paste(j, profile.list[i,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	                                names(prof.temp) <- c("temp.IID", "temp.SCORE","temp.PHENO")
	                                if(prev.files){
	                                    prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	                                    prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	                                    prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	                                }
	                                prof.out <- prof.temp
	                                names(prof.out) <- c("IID","SCORE","PHENO")
	                                prof <- prof.out
	                                prev.files <- T
	                            }
	                        }
	                    }
	                    if(geno.as.list){
	                        for(j in 1:22){
	                            if(paste(j, profile.list[i,1], sep = "") %in% list.files()){
	                                prof.temp <- read.table(paste(j, profile.list[i,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	                                names(prof.temp) <- c("temp.IID", "temp.SCORE","temp.PHENO")
	                                if(prev.files){
	                                    prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	                                    prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	                                    prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	                                }
	                                prof.out <- prof.temp
	                                names(prof.out) <- c("IID","SCORE","PHENO")
	                                prof <- prof.out
	                                prev.files <- T         
	                            }
	                            prev.files <- F
	                        }
	                    }
	                    if((!dosage & !geno.as.list)  | (dosage & is.na(dos.path.to.lists) & is.na(dos.list.file))){			
		                    ## Default - no fastscore
		                    prof <- read.table(paste(reduced.list[i], sep = ""), head  =T)
	                        prof <- prof[,c("IID", "PHENO", "SCORE")]
	                    }
	                    prof <- subset(prof, select=-c(PHENO))
	                    temp.pheno.data <- pheno.data[,c("ID", target.phenotypes[k])]
	                    names(temp.pheno.data) <- c("V1", "V2")
	                    prof <- merge(x = prof, by.x = "IID", y = temp.pheno.data, by.y = "V1")
	                    if(!no.regression){
	                        if(covary){
	                            prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
                                prof <- prof[prof$V2 != "" ,]
                                prof <- prof[prof$V2 != -9 ,]
                                prof <- prof[!is.na(prof$V2) ,]
	                            model.logit <- glm(V2 ~., family="binomial", data = prof[,c("V2","SCORE", covariates)])
	                            model.null <-  glm(V2 ~., family="binomial", data = prof[,c("V2", covariates)])
	                            p.out[i]  <- summary(model.logit)$coefficients[2,4]
	                            if(for.meta){
	                                coefficient[i]  <- summary(model.logit)$coefficients[2,1]
	                                s.err[i] <- summary(model.logit)$coefficients[2,2]
	                            }
		                        r2.out[i] <- NagelkerkeR2(model.logit)$R2 - NagelkerkeR2(model.null)$R2  
	                        }
	                        if(!covary){
                                prof <- prof[prof$V2 != "" ,]
                                prof <- prof[prof$V2 != -9 ,]
                                prof <- prof[!is.na(prof$V2) ,]
	                            p.out[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="binomial")))$coefficients[2,4]
	                            if(for.meta){
	                                coefficient[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="binomial")))$coefficients[2,1]
	                                s.err[i] <- summary(with(prof, glm(V2  ~ SCORE, family="binomial")))$coefficients[2,2]
	                            }
	                            r2.out[i] <- NagelkerkeR2(with(prof, glm(V2 ~ SCORE, family="binomial")))$R2
	                        }
	                    }
	                    prof.red <- prof[,c("IID","SCORE")]
	                    names(prof.red) <- c("IID", paste("pT", lists$V1[i], sep = "_"))
	                    if(report.individual.scores){
    	                    if(i == 1){
    	                        prof.all.scores <- prof.red
    	                    }
    	                    if(i > 1){
    	                        prof.all.scores <- merge(x = prof.all.scores, by.x = "IID", y = prof.red, by.y = "IID")    
    	                    }
	                    }
	                    nsnps[i] <- dim(ranges.list[ranges.list$V2 < lists$V1[i],])[1]
	                    if(i %in% (seq(1, length(lists$V1), 1)[seq(1, length(lists$V1), 1) %in% ceiling(length(lists$V1)*seq(0, 1,0.1))] )){
	                        cat(paste("Regression Models: " ,round(100*(i/length(lists$V1)),digits=-1), "% Complete \n",sep=""))
	                    }
                        if(calculate.abc.r2){
          	                abc.r2[i] <- estimateVg2FromP(prevalence1=prevalence1,
                                                            prevalence2=prevalence2,
                                                      	    n1=n1, 
                                                      	    sampling1= sampling1, 
                                                      	    n2= n2, 
                                                      	    sampling2= sampling2, 
                                                      	    nsnp= nsnps[i], 
                                                      	    plower = 0, 
                                                      	    pupper = lists$V1[i], 
                                                      	    binary=binary.target, 
                                                      	    p=p.out[i],
                                                      	    nullfraction=pi0
                                                      	    )$vg
                        }
	                }
	            }	
	            ## linear regression on phenotype
	            if(!binary.target){
	                for(i in 1:length(lists$V1)){
	                    if(dosage & (!is.na(dos.path.to.lists) | !is.na(dos.list.file))){
	                        for(j in 1:22){
	                            prof.temp <- read.table(paste(j, profile.list[j,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	                            names(prof.temp) <- c("temp.IID", "temp.SCORE")
	                            if(j > 1){
	                                prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	                                prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	                                prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	                            }
	                            prof.out <- prof.temp
	                            names(prof.out) <- c("IID","SCORE","PHENO")
	                        }
	                        prof <- prof.out
	                    }
	                    if(!dosage | (is.na(dos.path.to.lists) & is.na(dos.list.file))){			
	                        ## Default - no fastscore
	                        prof <- read.table(paste(reduced.list[i], sep = ""), head  =T)
	                        prof <- prof[,c("IID", "PHENO", "SCORE")]
	                    }
            	        prof <- subset(prof, select=-c(PHENO))
	                    temp.pheno.data <- pheno.data[,c("ID", target.phenotypes[k])]
	                    names(temp.pheno.data) <- c("V1", "V2")
	                    prof <- merge(x = prof, by.x = "IID", y = temp.pheno.data, by.y = "V1")
	                    if(!no.regression){
	                        if(covary){
	                            prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
                                prof <- prof[prof$V2 != "" ,]
                                prof <- prof[!is.na(prof$V2) ,]
		                        model.logit <- glm(V2  ~ ., family="gaussian", data = prof[,c("V2", "SCORE", covariates)])
		                        model.null <-  glm(V2  ~ ., family="gaussian", data = prof[,c("V2", covariates)])
		                        p.out[i]  <- summary(model.logit)$coefficients[2,4]
	                            if(for.meta){
	                                coefficient[i]  <- summary(model.logit)$coefficients[2,1]
	                                s.err[i] <- summary(model.logit)$coefficients[2,2]
	                            }
		                        r2.out[i] <- (1 - (sum( (prof$V2[complete.cases(prof[,c("V2", "SCORE", covariates)])]-predict(model.logit))^2 ) / sum( (prof$V2[complete.cases(prof[,c("V2", "SCORE", covariates)])]-mean(prof$V2[complete.cases(prof[,c("V2", "SCORE", covariates)])]))^2 ) ) )  - (1 - ( sum( (prof$V2[complete.cases(prof[,c("V2", "SCORE", covariates)])]-predict(model.null))^2 ) / sum( (prof$V2[complete.cases(prof[,c("V2", "SCORE", covariates)])]-mean(prof$V2[complete.cases(prof[,c("V2", "SCORE", covariates)])]))^2 ) ) )  
		                    }
	                        if(!covary){
                                prof <- prof[prof$V2 != "" ,]
                                prof <- prof[!is.na(prof$V2) ,]
	                            p.out[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="gaussian")))$coefficients[2,4]
	                            if(for.meta){
	                                coefficient[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="gaussian")))$coefficients[2,1]
	                                s.err[i] <- summary(with(prof, glm(V2  ~ SCORE, family="gaussian")))$coefficients[2,2]
	                            }
	                            r2.out[i] <-  1 - (sum( (prof$V2[complete.cases(prof[,c("V2", "SCORE")])]-predict(with(prof[,c("V2", "SCORE")], glm(V2 ~ SCORE, family="gaussian"))))^2 )/sum( (prof$V2[complete.cases(prof[,c("V2", "SCORE")])]-mean(prof$V2[complete.cases(prof[,c("V2", "SCORE")])]))^2))   			
	                        }
	                    }
	                    prof.red <- prof[,c("IID","SCORE")]
	                    names(prof.red) <- c("IID", paste("pT", lists$V1[i], sep = "_"))
	                    if(report.individual.scores){
	                        if(i == 1){
	                            prof.all.scores <- prof.red
	                        }
	                        if(i > 1){
	                            prof.all.scores <- merge(x = prof.all.scores, by.x = "IID", y = prof.red, by.y = "IID")    
	                        }
	                    }
	                    nsnps[i] <- dim(ranges.list[ranges.list$V2 < lists$V1[i],])[1]
	                    if(i %in% (seq(1, length(lists$V1), 1)[seq(1, length(lists$V1), 1) %in% ceiling(length(lists$V1)*seq(0, 1,0.1))] )){
	                        cat(paste("Regression Models: " ,round(100*(i/length(lists$V1)),digits=-1), "% Complete \n",sep=""))
	                    }
                        if(calculate.abc.r2){
          	                abc.r2[i] <- estimateVg2FromP(prevalence1=prevalence1,
                                                  	        prevalence2=prevalence2,
                                                  	        n1=n1, 
                                                  	        sampling1= sampling1, 
                                                  	        n2= n2, 
                                                  	        sampling2= sampling2, 
                                                  	        nsnp= nsnps[i], 
                                                  	        plower = 0, 
                                                  	        pupper = lists$V1[i], 
                                                  	        binary=binary.target, 
                                                  	        p=p.out[i],
                                                  	        nullfraction=pi0
                                                  	        )$vg
                        }
	                }
	            }
	            top.thresh[k] <- lists$V1[which.min(p.out)]
  	            top.thresh.pval[k] <- p.out[which.min(p.out)]
  	            top.thresh.r2[k] <- r2.out[which.min(p.out)]
  	    
  	    
	            thresh <- lists$V1 
                if(!calculate.abc.r2){
	                if(!for.meta){
	                    output <- data.frame(thresh, p.out, r2.out, nsnps)
	                }
	                if(for.meta){
	                    output <- data.frame(thresh, p.out, r2.out, nsnps, coefficient, s.err)
	                }
	            }
                if(calculate.abc.r2){
	                if(!for.meta){
	                    output <- data.frame(thresh, p.out, abc.r2, nsnps)
	                }
	                if(for.meta){
	                    output <- data.frame(thresh, p.out,abc.r2, nsnps, coefficient, s.err)
	                }
	            }
	            if(!multiple.base.phenotypes){
  	                write.table(output, paste(figname, target.phenotypes[k], "RAW_RESULTS_DATA.txt", sep = "_"), col.names = T, row.names = F, quote  =F)
	            }
	            if(multiple.base.phenotypes){
  	                write.table(output, paste(figname, base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k], "RAW_RESULTS_DATA.txt", sep = "_"), col.names = T, row.names = F, quote  =F)
	            }
	            if(report.individual.scores){
	                if(!report.best.score.only & !multiple.base.phenotypes){
	                    write.table(prof.all.scores, paste(figname, target.phenotypes[k], "SCORES_AT_ALL_THRESHOLDS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	                }
	                if(!report.best.score.only & multiple.base.phenotypes){
	                    write.table(prof.all.scores, paste(figname, base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],   "SCORES_AT_ALL_THRESHOLDS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	                }
	                if(report.best.score.only & !multiple.base.phenotypes){
	                    write.table(prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))], paste(figname, target.phenotypes[k], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	                }
	                if(report.best.score.only & multiple.base.phenotypes){
	                    write.table(prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))], paste(figname,base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	                }
	            }
	        }
	        ## end of loop through target phenotypes
	        if(!multiple.base.phenotypes){  
    	        write.table(data.frame(target.phenotypes, target.phenotypes.binary, top.thresh, top.thresh.pval, top.thresh.r2), paste(figname, "TOP_SCORES_ACROSS_PHENOTYPES.txt",sep="_"), col.names=T, row.names=F, quote=F)
            }	
	        if(multiple.base.phenotypes){  
                write.table(data.frame(target.phenotypes, target.phenotypes.binary, top.thresh, top.thresh.pval, top.thresh.r2), paste(figname, base.phenotypes.names[basePhen], "PREDICTING_TOP_SCORES_ACROSS_TARGET_PHENOTYPES.txt",sep="_"), col.names=T, row.names=F, quote=F)
            }	

        }
  

	    if(no.regression & !multiple.target.phenotypes){
	        if(report.individual.scores & !multiple.base.phenotypes){
	            write.table(prof.all.scores, paste(figname, "SCORES_AT_ALL_THRESHOLDS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	        }
	        if(report.individual.scores & multiple.base.phenotypes){
	            write.table(prof.all.scores, paste(figname,base.phenotypes.names[basePhen], "SCORES_AT_ALL_THRESHOLDS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	        }
	        if(cleanup){
	            cat(" ################################# \n # \n #   Cleanup \n # \n ################################# \n")
	            system("rm PROFILES.0.5*.profile", ignore.stdout=T,ignore.stderr=T)	
	            system("rm PROFILES.0.4*.profile", ignore.stdout=T,ignore.stderr=T)
	            system("rm PROFILES.0.3*.profile", ignore.stdout=T,ignore.stderr=T)
	            system("rm PROFILES.0.2*.profile", ignore.stdout=T,ignore.stderr=T)
	            system("rm PROFILES.0.1*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm PROFILES.0.0*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm PROFILES.*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm PROFILES.S1*.profile", ignore.stdout=T,ignore.stderr=T)	
        	    system("rm PROFILES.S2*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm PROFILES.S3*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm PROFILES.S4*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm PROFILES.S5*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm PROFILES.S6*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm PROFILES.S7*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm PROFILES.S8*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm PROFILES.S9*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.S*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.0.5*.profile", ignore.stdout=T,ignore.stderr=T)	
        	    system("rm *PROFILES.0.4*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.0.3*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.0.2*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.0.1*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.0.0*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.S1*.profile", ignore.stdout=T,ignore.stderr=T)	
        	    system("rm *PROFILES.S2*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.S3*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.S4*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.S5*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.S6*.profile", ignore.stdout=T,ignore.stderr=T) 
        	    system("rm *PROFILES.S7*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.S8*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.S9*.profile", ignore.stdout=T,ignore.stderr=T)
        	    system("rm *PROFILES.S*.profile", ignore.stdout=T,ignore.stderr=T) 
        	    system("rm flip*", ignore.stdout=T,ignore.stderr=T)
        	    system("rm clean*", ignore.stdout=T,ignore.stderr=T)
        	    system("rm profile_list", ignore.stdout=T,ignore.stderr=T)
        	    system("rm rangelist.txt", ignore.stdout=T,ignore.stderr=T)
        	    system("rm rangelist_ranges", ignore.stdout=T,ignore.stderr=T)
        	    system("rm rawfile.raw", ignore.stdout=T,ignore.stderr=T)
        	    system("rm reordered_base", ignore.stdout=T,ignore.stderr=T)
        	    system("rm target_no_mhc*", ignore.stdout=T,ignore.stderr=T)
        	    system("rm ANCESTRY_INFORMATIVE_DIMENSIONS.cluster*", ignore.stdout=T,ignore.stderr=T)
        	    system("rm reordered_base", ignore.stdout=T,ignore.stderr=T)
        	    system("rm ANCESTRY_INFORMATIVE_DIMENSIONS.log", ignore.stdout=T,ignore.stderr=T)
        	    system("rm head_disc", ignore.stdout=T,ignore.stderr=T)
        	    system("rm LE_SNPs", ignore.stdout=T,ignore.stderr=T)
        	    system("rm prune_target*", ignore.stdout=T,ignore.stderr=T)
        	    system("rm plink.log", ignore.stdout=T,ignore.stderr=T)
        	    system("rm HEADER", ignore.stdout=T,ignore.stderr=T)
        	    system("rm PROFILES.hh", ignore.stdout=T,ignore.stderr=T)
        	    system("rm PROFILES.log", ignore.stdout=T,ignore.stderr=T)
        	    system("rm synonymous_snps*", ignore.stdout=T,ignore.stderr=T)
        	    system("rm synonymous_snps", ignore.stdout=T,ignore.stderr=T) 
        	    system("rm non_synonymous_snps_only*", ignore.stdout=T,ignore.stderr=T)
        	    system("rm mhc.txt", ignore.stdout=T,ignore.stderr=T)
        	    system("rm Rplots.pdf", ignore.stdout=T,ignore.stderr=T)
        	    system("rm temp.raw", ignore.stdout=T,ignore.stderr=T)
        	    system("rm TARGET_SNPs", ignore.stdout=T,ignore.stderr=T)
        	    system("rm base_SNPS", ignore.stdout=T,ignore.stderr=T)  
        	    system("rm Complete_Allele_List.txt", ignore.stdout=T,ignore.stderr=T)
	        }
	        if(print.time){
	            cat(" ################################# \n # \n #   Print time \n # \n ################################# \n")
	        }
	        running.time <- proc.time()[3] - start.time 
	        if(running.time < 60){
	            out.run.time <- round(running.time, digits = 2)
	            out.time.units <- "seconds"
	        }
	        if(running.time > 60 & running.time < 3600){
	            out.run.time <- round((running.time/60), digits = 2)
	            out.time.units <- "minutes"
	        }
	        if(running.time > 3600){
	            out.run.time <- round((running.time/3600), digits = 2)
	            out.time.units <- "hours"
	        }
	        if(print.time){
	            print(paste("RUNNING TIME: ", 	out.run.time, out.time.units, sep = " "))
	        }
	        quit()
        }
	
	    if(!multiple.target.phenotypes){
	        thresh <- lists$V1 
            if(!calculate.abc.r2){
	            if(!for.meta){
	                output <- data.frame(thresh, p.out, r2.out, nsnps)
	            }
	            if(for.meta){
	                output <- data.frame(thresh, p.out, r2.out, nsnps, coefficient, s.err)
	            }
	        }
            if(calculate.abc.r2){
	            if(!for.meta){
	                output <- data.frame(thresh, p.out, abc.r2, nsnps)
	            }
	            if(for.meta){
	                output <- data.frame(thresh, p.out,abc.r2, nsnps, coefficient, s.err)
	            }
	        }
	        if(!multiple.base.phenotypes){
  	            write.table(output, paste(figname, "RAW_RESULTS_DATA.txt", sep = "_"), col.names = T, row.names = F, quote  =F)
  	        }
	        if(multiple.base.phenotypes){
  	            write.table(output, paste(figname,base.phenotypes.names[basePhen], "RAW_RESULTS_DATA.txt", sep = "_"), col.names = T, row.names = F, quote  =F)
  	        }
  	  
	        if(report.individual.scores){
	            if(!report.best.score.only & !multiple.base.phenotypes){
	                write.table(prof.all.scores, paste(figname, "SCORES_AT_ALL_THRESHOLDS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	            }
	            if(!report.best.score.only & multiple.base.phenotypes){
	                write.table(prof.all.scores, paste(figname,base.phenotypes.names[basePhen],"SCORES_AT_ALL_THRESHOLDS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	            }
	            if(report.best.score.only & !multiple.base.phenotypes){
	                write.table(prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))], paste(figname, "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	            }
	            if(report.best.score.only & multiple.base.phenotypes){
	                write.table(prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))], paste(figname,base.phenotypes.names[basePhen], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	            }
	        }
	    }
    
	    if(!multiple.base.phenotypes){
	        if(!multiple.target.phenotypes){
                if(quantiles & ggfig){
                    cat(" ################################# \n # \n #   Quantiles Plots \n # \n ################################# \n")
                    if(binary.target){
                        if(report.individual.scores & !report.best.score.only){
                            scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
                        }
                        if(report.individual.scores & report.best.score.only){
                            scores.internal <- read.table(paste(figname, "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
                        }
                        names(scores.internal) <- c("ID", "SCORE")
                        for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
                        for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
                        for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))

                        if(covary == F){
                            or.quantiles <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1])
                            ci.quantiles.u <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
                            ci.quantiles.l <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
                        }
                        if(covary == T){
                            for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")
                            or.quantiles <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1])
                            ci.quantiles.u <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
                            ci.quantiles.l <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
                        }
                        or.quantiles[1] <- 1 
                        ci.quantiles.u[1] <- 1
                        ci.quantiles.l[1] <- 1
                        quant.list <- seq(1, num.quantiles, 1)
                        quant.list <- quant.list[quant.list != quant.ref]     
                        quantiles.for.table <- c(quant.ref, quant.list)
                        quantiles.df <- data.frame(or.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
                        names(quantiles.df) <- c("OR", "CI.U", "CI.L", "DEC")
                        quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
                        quantiles.plot <- ggplot(quantiles.df) + 
                        geom_point(aes(x = DEC, y = OR), colour = "royalblue2", size=4) + 
                        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                        axis.line = element_line(colour = "black",size=0.5)) + 
                        geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = OR, x = DEC), colour = "royalblue2", size = 0.9) + 
                        ylab("Odds Ratio for Score on Phenotype") + 
                        xlab("Quantiles for Polygenic Score") + 
                        scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
                    }  
                    if(!binary.target){
                        if(report.individual.scores & !report.best.score.only){
                            scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
                        }
                        if(report.individual.scores & report.best.score.only){
                            scores.internal <- read.table(paste(figname, "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
                        }
                        names(scores.internal) <- c("ID", "SCORE")
                        for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
                        for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
                        for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))
    
                        if(covary == F){
                            coef.quantiles <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1]
                            ci.quantiles.u <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
                            ci.quantiles.l <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
                        }
                        if(covary == T){
                            for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")
                            coef.quantiles <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1]
                            ci.quantiles.u <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
                            ci.quantiles.l <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
                        }
                        coef.quantiles[1] <- 0 
                        ci.quantiles.u[1] <- 0
                        ci.quantiles.l[1] <- 0
                        quant.list <- seq(1, num.quantiles, 1)
                        quant.list <- quant.list[quant.list != quant.ref]     
                        quantiles.for.table <- c(quant.ref, quant.list)
                        quantiles.df <- data.frame(coef.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
                        names(quantiles.df) <- c("Coef", "CI.U", "CI.L", "DEC")
                        quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
                        quantiles.plot <- ggplot(quantiles.df) + 
                          geom_point(aes(x = DEC, y = Coef), colour = "royalblue2", size=4) + 
                          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                          axis.line = element_line(colour = "black",size=0.5)) + 
                          geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = Coef, x = DEC), colour = "royalblue2", size = 0.9) + 
                          ylab("Change in Phenotype given score in quantiles") + 
                          xlab("Quantiles for Polygenic Score") + 
                          scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
                    }  
                    quantiles.plot
                    ggsave(paste(figname, "QUANTILES_PLOT.png", sep = "_"))  
                }	

        	    if(!fastscore){
        	        cat(" ################################# \n # \n #   High Density Plots \n # \n ################################# \n")
        	    }
	            
        	    if(best.thresh.on.bar){
        	        #add max thresh
        	        barchart.levels <- c(barchart.levels.old, output$thresh[which.min(output$p.out)])
        	        barchart.levels <- barchart.levels[!duplicated(barchart.levels)]	
        	        barchart.levels <- sort(barchart.levels, decreasing = F)
        	    }
	            if(!fastscore & !ggfig){
	                if(!scatter.R2){
	                    png(paste(figname,"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	                    with(output, plot(	x = output$thresh,
	                                        y = -log10(output$p.out),
                            	          xlab = expression(italic(P)-value~threshold),
                            	          ylab = bquote(-log[10]~model~italic(P)-value),
                            	          pch = 19,
                            	          type = "b"
	                                    ))
	                    with(output[output$thresh %in% barchart.levels.old,], points(thresh, -log10(p.out), col = "green", pch = 19, type= "b"))
	                }
        	        if(scatter.R2){
        	            png(paste(figname,"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
        	            with(output, plot(	x = output$thresh,
        	                y = r2.out,
        	                xlab = expression(italic(P)-value~threshold),
        	                ylab = bquote(-log[10]~model~italic(P)-value),
        	                pch = 19,
        	                type = "b"
        	            ))
        	            with(output[output$thresh %in% barchart.levels.old,], points(thresh, r2.out, col = "green", pch = 19, type= "b"))
        	        }	
	                dev.off()
	            }
	            if(!fastscore & ggfig){
	                if(!scatter.R2){
            	        ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = -log10(p.out))) +
            	          geom_line(aes(x = thresh, y = -log10(p.out))) +
                          xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
                          ylab(bquote(PRS~model~fit:~italic(P)-value~(-log[10]))) +
            	          geom_line(aes(thresh,  -log10(p.out)), colour = "green",
            	          data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
            	          geom_point(aes(thresh,  -log10(p.out)), colour = "green",
            	         data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
            	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                      axis.line = element_line(colour = "black",size=0.5))
            	        ggfig.points
            	        ggsave(paste(figname,"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	                }
	                if(scatter.R2){
            	        ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = r2.out)) +
            	          geom_line(aes(x = thresh, y = r2.out)) +
                          xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
            	          geom_line(aes(thresh,  r2.out), colour = "green",
            	          data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
            	          geom_point(aes(thresh,  r2.out), colour = "green",
            	          data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
            	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                      axis.line = element_line(colour = "black",size=0.5))
                        if(binary.target){
                       	    if(!calculate.abc.r2){
                                ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
                            }
                       	    if(calculate.abc.r2){
                                ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
                            }
                        }
                        if(!binary.target){
                            ggfig.points <- ggfig.points + ylab(paste("PRS model fit:  ",expression(paste(R^2))))
                        }
            	        ggfig.points
            	        ggsave(paste(figname,"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	                }	 
	            }
	            options(scipen=0,digits=7)
	            output <- read.table(paste(figname, "RAW_RESULTS_DATA.txt", sep = "_"),head=T)
	            output <- output[,c(1:3)]
	            names(output) <- c("thresh", "p.out", "r2.out")
	            output <- output[output$thresh %in% barchart.levels,]
	            output$print.p[round(output$p.out, digits = 3) != 0] <- round(output$p.out[round(output$p.out, digits = 3) != 0], digits = 3)
	            output$print.p[round(output$p.out, digits = 3) == 0] <- format(output$p.out[round(output$p.out, digits = 3) == 0], digits=2)
	            output$print.p <- sub("e", "*x*10^", output$print.p)
	            cat(" ################################# \n # \n #   Barplots \n # Bars for inclusion can be changed using the barchart.levels option \n ################################# \n")
	            if(ggfig){
	                if(!bar.col.is.pval){
	                    ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = factor(thresh)), stat="identity") +     scale_fill_brewer(palette=barpalatte, name = expression(italic(P)-value~threshold)) + xlab(expression(italic(P)-value~threshold~(italic(P)[T])))   
	                }
	                if(bar.col.is.pval){
	                    ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = -log10(p.out)), stat="identity") +     scale_fill_gradient(low= bar.col.is.pval.lowcol, high= bar.col.is.pval.highcol, name =bquote(atop(-log[10]~model,italic(P)-value),)) +  xlab(expression(italic(P)-value~threshold~(italic(P)[T]))) 	      
	                }
	                ggfig.plot <- ggfig.plot + geom_text(aes(x = factor(thresh), y = r2.out, label = paste(print.p)), vjust = -1.5, hjust = 0, angle = 45, cex = 2.8, parse=T) +
	                                scale_y_continuous(limits = c(0, max(output$r2.out)*1.25)) +
	                                theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                    axis.line = element_line(colour = "black",size=0.5))
                    if(binary.target){
   	                    if(!calculate.abc.r2){
                            ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
                        }
   	                    if(calculate.abc.r2){
                            ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
                        }
                    }
                    if(!binary.target){
                        ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ",R^2)))
                    }
	                ggfig.plot
	                ggsave(paste(figname, "_BARPLOT_", Sys.Date(), ".png", sep = ""))		
	            }  
	            if(!ggfig) {
	                png(paste(figname, "_BARPLOT_", Sys.Date(), ".png", sep = ""))
	                plot.fig  <- with(output, barplot(r2.out, 
                    	        names = thresh,
                    	        main = "", 
                    	        col = "red",
                    	        xlab = expression(italic(P)[T]),
                    	        ylab = expression(R^2),
                    	        ylim = c(0,  max(output$r2.out)*1.25)  ))
                    	        text( parse(text=paste(
                    	        output$print.p)), 
                    	        x = plot.fig+0.1, 
                    	        y =  output$r2.out+ (max(output$r2.out)*1.125-max(output$r2.out)), 
                    	        srt = 45)
        	        dev.off()
	            }
	        }
	        if(multiple.target.phenotypes){
	            for(k in 1:length(target.phenotypes)){
	                output <- read.table(paste(figname, target.phenotypes[k], "RAW_RESULTS_DATA.txt", sep = "_"), head=T)
	      
	                if(quantiles & ggfig){
                        cat(" ################################# \n # \n #   Quantiles Plots \n # \n ################################# \n")
                        if(binary.target){
                            if(report.individual.scores & !report.best.score.only){
                                scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
                            }
                            if(report.individual.scores & report.best.score.only){
                                scores.internal <- read.table(paste(figname,target.phenotypes[k], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
                            }
                            names(scores.internal) <- c("ID", "SCORE")
                            for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
                            for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
                            for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))

                            if(covary == F){
                                or.quantiles <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1])
                                ci.quantiles.u <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
                                ci.quantiles.l <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
                            }
                            if(covary == T){
                                for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")

                                or.quantiles <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1])
                                ci.quantiles.u <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2]))
                                ci.quantiles.l <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2]))
                            }
                            or.quantiles[1] <- 1 
                            ci.quantiles.u[1] <- 1
                            ci.quantiles.l[1] <- 1
                            quant.list <- seq(1, num.quantiles, 1)
                            quant.list <- quant.list[quant.list != quant.ref]     
                            quantiles.for.table <- c(quant.ref, quant.list)
                            quantiles.df <- data.frame(or.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
                            names(quantiles.df) <- c("OR", "CI.U", "CI.L", "DEC")
                            quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
                            quantiles.plot <- ggplot(quantiles.df) + 
                              geom_point(aes(x = DEC, y = OR), colour = "royalblue2", size=4) + 
                              theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                              axis.line = element_line(colour = "black",size=0.5)) + 
                              geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = OR, x = DEC), colour = "royalblue2", size = 0.9) + 
                              ylab("Odds Ratio for Score on Phenotype") + 
                              xlab("quantiles for Polygenic Score") + 
                              scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
                        }  
                        if(!binary.target){
                            if(report.individual.scores & !report.best.score.only){
                                scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
                            }
                            if(report.individual.scores & report.best.score.only){
                                scores.internal <- read.table(paste(figname,target.phenotypes[k], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
                            }
                            names(scores.internal) <- c("ID", "SCORE")
                            for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
                            for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
                            for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))

                            if(covary == F){
                                coef.quantiles <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1]
                                ci.quantiles.u <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
                                ci.quantiles.l <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
                            }
                            if(covary == T){
                                for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")
                                coef.quantiles <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1]
                                ci.quantiles.u <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
                                ci.quantiles.l <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
                            }
                            coef.quantiles[1] <- 0 
                            ci.quantiles.u[1] <- 0
                            ci.quantiles.l[1] <- 0
                            quant.list <- seq(1, num.quantiles, 1)
                            quant.list <- quant.list[quant.list != quant.ref]     
                            quantiles.for.table <- c(quant.ref, quant.list)
                            quantiles.df <- data.frame(coef.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
                            names(quantiles.df) <- c("Coef", "CI.U", "CI.L", "DEC")
                            quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
                            quantiles.plot <- ggplot(quantiles.df) + 
                              geom_point(aes(x = DEC, y = Coef), colour = "royalblue2", size=4) + 
                              theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                              axis.line = element_line(colour = "black",size=0.5)) + 
                              geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = Coef, x = DEC), colour = "royalblue2", size = 0.9) + 
                              ylab("Change in Phenotype given score in quantiles") + 
                              xlab("Quantiles for Polygenic Score") + 
                              scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
                        }  
                        quantiles.plot
                        ggsave(paste(figname,target.phenotypes[k], "QUANTILES_PLOT.png", sep = "_"))  
                    }	

	      
        	        if(!fastscore){
        	            cat(" ################################# \n # \n #   High Density Plots \n # \n ################################# \n")
        	        }
        	        if(best.thresh.on.bar){
        	            #add max thresh
        	            barchart.levels <- c(barchart.levels.old, output$thresh[which.min(output$p.out)])
        	            barchart.levels <- barchart.levels[!duplicated(barchart.levels)]	
        	            barchart.levels <- sort(barchart.levels, decreasing = F)
        	        }
	                if(!fastscore & !ggfig){
            	        if(!scatter.R2){
                	        png(paste(figname,"_",target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
                	        with(output, plot(	x = output$thresh,
                	            y = -log10(output$p.out),
                	            xlab = expression(italic(P)-value~threshold),
                	            ylab = bquote(-log[10]~model~italic(P)-value),
                	            pch = 19,
                	            type = "b"
                	        ))
                	        with(output[output$thresh %in% barchart.levels.old,], points(thresh, -log10(p.out), col = "green", pch = 19, type= "b"))
                        }
        	            if(scatter.R2){
        	                png(paste(figname,"_",target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
        	                with(output, plot(	x = output$thresh,
                	            y = r2.out,
                	            xlab = expression(italic(P)-value~threshold),
                	            ylab = bquote(-log[10]~model~italic(P)-value),
                	            pch = 19,
                	            type = "b"
        	                ))
        	                with(output[output$thresh %in% barchart.levels.old,], points(thresh, r2.out, col = "green", pch = 19, type= "b"))
        	            }	
        	            dev.off()
	                }  
        	        if(!fastscore & ggfig){
        	            if(!scatter.R2){
        	                ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = -log10(p.out))) +
                	            geom_line(aes(x = thresh, y = -log10(p.out))) +
                                xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
                                ylab(bquote(PRS~model~fit:~italic(P)-value~(-log[10]))) +
                	            geom_line(aes(thresh,  -log10(p.out)), colour = "green",
                	             data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
                	             geom_point(aes(thresh,  -log10(p.out)), colour = "green",
                	            data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
                	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                          axis.line = element_line(colour = "black",size=0.5))
                	          ggfig.points
            	            ggsave(paste(figname,"_",target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
        	            }
	                    if(scatter.R2){
            	            ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = r2.out)) +
                	            geom_line(aes(x = thresh, y = r2.out)) +
                                xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
                	            geom_line(aes(thresh,  r2.out), colour = "green",
                	            data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
                	            geom_point(aes(thresh,  r2.out), colour = "green",
                	            data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
                	            theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black",size=0.5))
                            if(binary.target){
                           	    if(!calculate.abc.r2){
                                    ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
                                }
                           	    if(calculate.abc.r2){
                                    ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
                                }
                            }
                            if(!binary.target){
                                ggfig.points <- ggfig.points + ylab(paste("PRS model fit:  ",expression(paste(R^2))))
                            }
            	            ggfig.points
            	            ggsave(paste(figname,"_",target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	                    }	 
	                }
	                options(scipen=0,digits=7)
	                output <- output[,c(1:3)]
	                names(output) <- c("thresh", "p.out", "r2.out")
	                output <- output[output$thresh %in% barchart.levels,]
            	    output$print.p[round(output$p.out, digits = 3) != 0] <- round(output$p.out[round(output$p.out, digits = 3) != 0], digits = 3)
            	    output$print.p[round(output$p.out, digits = 3) == 0] <- format(output$p.out[round(output$p.out, digits = 3) == 0], digits=2)
            	    output$print.p <- sub("e", "*x*10^", output$print.p)
            	    cat(" ################################# \n # \n #   Barplots \n # Bars for inclusion can be changed using the barchart.levels option \n ################################# \n")
	                if(ggfig){
        	            if(!bar.col.is.pval){
        	                ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = factor(thresh)), stat="identity") +     scale_fill_brewer(palette=barpalatte, name = expression(italic(P)-value~threshold)) + theme(axis.text.x=element_blank(), axis.title.x=element_blank())  +   xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  
        	            }
	                    if(bar.col.is.pval){
            	            ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = -log10(p.out)), stat="identity") +     scale_fill_gradient(low= bar.col.is.pval.lowcol, high= bar.col.is.pval.highcol, name =bquote(atop(-log[10]~model,italic(P)-value),)) +  xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  
	                    }
	                        ggfig.plot <- ggfig.plot + geom_text(aes(x = factor(thresh), y = r2.out, label = paste(print.p)), vjust = -1.5, hjust = 0, angle = 45, cex = 2.8, parse=T) +
	                                        scale_y_continuous(limits = c(0, max(output$r2.out)*1.25))  + 
	                                        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                            axis.line = element_line(colour = "black",size=0.5))
                        if(binary.target){
                   	        if(!calculate.abc.r2){
                                ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
                            }
                   	        if(calculate.abc.r2){
                                ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
                            }
                        }
                        if(!binary.target){
                            ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ",R^2)))
                        }
        	            ggfig.plot
        	            ggsave(paste(figname,"_",target.phenotypes[k],"_BARPLOT_", Sys.Date(), ".png", sep = ""))		
	                }  
        	        if(!ggfig) {
            	        png(paste(figname,"_",target.phenotypes[k],"_BARPLOT_", Sys.Date(), ".png", sep = ""))
            	        plot.fig  <- with(output, barplot(r2.out, 
            	          names = thresh,
            	          main = "", 
            	          col = "red",
            	          xlab = expression(italic(P)[T]),
            	          ylab = expression(R^2),
            	          ylim = c(0,  max(output$r2.out)*1.25)  ))
            	          text( parse(text=paste(
            	          output$print.p)), 
            	          x = plot.fig+0.1, 
            	          y =  output$r2.out+ (max(output$r2.out)*1.125-max(output$r2.out)), 
            	          srt = 45)
            	        dev.off()
        	        }
	            }
	        }
	    }
	    if(multiple.base.phenotypes){
    	    if(!multiple.target.phenotypes){
    	        output <- read.table(paste(figname,base.phenotypes.names[basePhen], "RAW_RESULTS_DATA.txt", sep = "_"), head=T)
                if(quantiles & ggfig){
                    cat(" ################################# \n # \n #   Quantiles Plots \n # \n ################################# \n")
                    if(binary.target){
                        if(report.individual.scores & !report.best.score.only){
                            scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
                        }
                        if(report.individual.scores & report.best.score.only){
                            scores.internal <- read.table(paste(figname,base.phenotypes.names[basePhen], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
                        }
                        names(scores.internal) <- c("ID", "SCORE")
                        for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
                        for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
                        for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))
    
                        if(covary == F){
                            or.quantiles <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1])
                            ci.quantiles.u <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
                            ci.quantiles.l <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
                        }
                        if(covary == T){
                            for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")
    
                            or.quantiles <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1])
                            ci.quantiles.u <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2]))
                            ci.quantiles.l <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2]))
                        }
                        or.quantiles[1] <- 1 
                        ci.quantiles.u[1] <- 1
                        ci.quantiles.l[1] <- 1
                        quant.list <- seq(1, num.quantiles, 1)
                        quant.list <- quant.list[quant.list != quant.ref]     
                        quantiles.for.table <- c(quant.ref, quant.list)
                        quantiles.df <- data.frame(or.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
                        names(quantiles.df) <- c("OR", "CI.U", "CI.L", "DEC")
                        quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
                        quantiles.plot <- ggplot(quantiles.df) + 
                          geom_point(aes(x = DEC, y = OR), colour = "royalblue2", size=4) + 
                          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                          axis.line = element_line(colour = "black",size=0.5)) + 
                          geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = OR, x = DEC), colour = "royalblue2", size = 0.9) + 
                          ylab("Odds Ratio for Score on Phenotype") + 
                          xlab("quantiles for Polygenic Score") + 
                          scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
                    }  
                    if(!binary.target){
                        if(report.individual.scores & !report.best.score.only){
                            scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
                        }
                        if(report.individual.scores & report.best.score.only){
                            scores.internal <- read.table(paste(figname, base.phenotypes.names[basePhen],"SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
                        }
                        names(scores.internal) <- c("ID", "SCORE")
                        for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
                        for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
                        for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))
        
                        if(covary == F){
                            coef.quantiles <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1]
                            ci.quantiles.u <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
                            ci.quantiles.l <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
                        }
                        if(covary == T){
                            for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")
                            coef.quantiles <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1]
                            ci.quantiles.u <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
                            ci.quantiles.l <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
                        }
                        coef.quantiles[1] <- 0 
                        ci.quantiles.u[1] <- 0
                        ci.quantiles.l[1] <- 0
                        quant.list <- seq(1, num.quantiles, 1)
                        quant.list <- quant.list[quant.list != quant.ref]     
                        quantiles.for.table <- c(quant.ref, quant.list)
                        quantiles.df <- data.frame(coef.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
                        names(quantiles.df) <- c("Coef", "CI.U", "CI.L", "DEC")
                        quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
                        quantiles.plot <- ggplot(quantiles.df) + 
                             geom_point(aes(x = DEC, y = Coef), colour = "royalblue2", size=4) + 
                             theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                             axis.line = element_line(colour = "black",size=0.5)) + 
                             geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = Coef, x = DEC), colour = "royalblue2", size = 0.9) + 
                             ylab("Change in Phenotype given score in quantiles") + 
                             xlab("quantiles for Polygenic Score") + 
                             scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
                    }  
                    quantiles.plot
                    ggsave(paste(figname, base.phenotypes.names[basePhen],"QUANTILES_PLOT.png", sep = "_"))  
                }	
    
        	    if(!fastscore){
        	        cat(" ################################# \n # \n #   High Density Plots \n # \n ################################# \n")
        	    }
        	    if(best.thresh.on.bar){
        	        #add max thresh
        	        barchart.levels <- c(barchart.levels.old, output$thresh[which.min(output$p.out)])
        	        barchart.levels <- barchart.levels[!duplicated(barchart.levels)]	
        	        barchart.levels <- sort(barchart.levels, decreasing = F)
        	    }
    	        if(!fastscore & !ggfig){
    	            if(!scatter.R2){
    	                png(paste(figname,"_",base.phenotypes.names[basePhen],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
        	            with(output, plot(	x = output$thresh,
            	          y = -log10(output$p.out),
            	          xlab = expression(italic(P)-value~threshold),
            	          ylab = bquote(-log[10]~model~italic(P)-value),
            	          pch = 19,
            	          type = "b"
    	                ))
    	                with(output[output$thresh %in% barchart.levels.old,], points(thresh, -log10(p.out), col = "green", pch = 19, type= "b"))
    	            }
    	            if(scatter.R2){
            	        png(paste(figname,"_",base.phenotypes.names[basePhen],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
            	        with(output, plot(	x = output$thresh,
            	          y = r2.out,
            	          xlab = expression(italic(P)-value~threshold),
            	          ylab = bquote(-log[10]~model~italic(P)-value),
            	          pch = 19,
            	          type = "b"
            	        ))
    	                with(output[output$thresh %in% barchart.levels.old,], points(thresh, r2.out, col = "green", pch = 19, type= "b"))
    	            }	
    	            dev.off()
    	        }
    	        if(!fastscore & ggfig){
    	            if(!scatter.R2){
    	                ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = -log10(p.out))) +
            	          geom_line(aes(x = thresh, y = -log10(p.out))) +
                          xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
                          ylab(bquote(PRS~model~fit:~italic(P)-value~(-log[10]))) +
            	          geom_line(aes(thresh,  -log10(p.out)), colour = "green",
            	          data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
            	          geom_point(aes(thresh,  -log10(p.out)), colour = "green",
            	          data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
            	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                      axis.line = element_line(colour = "black",size=0.5))
            	        ggfig.points
        	            ggsave(paste(figname,"_",base.phenotypes.names[basePhen],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
    	            }
    	            if(scatter.R2){
            	        ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = r2.out)) +
            	          geom_line(aes(x = thresh, y = r2.out)) +
                          xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
            	          geom_line(aes(thresh,  r2.out), colour = "green",
            	          data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
            	          geom_point(aes(thresh,  r2.out), colour = "green",
            	         data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
            	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                        axis.line = element_line(colour = "black",size=0.5))
                        if(binary.target){
                       	    if(!calculate.abc.r2){
                                ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
                            }
           	                if(calculate.abc.r2){
                                ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
                            }
                        }
                        if(!binary.target){
                            ggfig.points <- ggfig.points + ylab(paste("PRS model fit:  ",expression(paste(R^2))))
                        }
        	            ggfig.points
        	            ggsave(paste(figname,"_",base.phenotypes.names[basePhen],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
    	            }	 
    	        }
    	        options(scipen=0,digits=7)
        	    output <- read.table(paste(figname,base.phenotypes.names[basePhen], "RAW_RESULTS_DATA.txt", sep = "_"),head=T)
        	    output <- output[,c(1:3)]
        	    names(output) <- c("thresh", "p.out", "r2.out")
        	    output <- output[output$thresh %in% barchart.levels,]
        	    output$print.p[round(output$p.out, digits = 3) != 0] <- round(output$p.out[round(output$p.out, digits = 3) != 0], digits = 3)
        	    output$print.p[round(output$p.out, digits = 3) == 0] <- format(output$p.out[round(output$p.out, digits = 3) == 0], digits=2)
        	    output$print.p <- sub("e", "*x*10^", output$print.p)
        	    cat(" ################################# \n # \n #   Barplots \n # Bars for inclusion can be changed using the barchart.levels option \n ################################# \n")
    	        if(ggfig){
    	            if(!bar.col.is.pval){
    	                ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = factor(thresh)), stat="identity") +     scale_fill_brewer(palette=barpalatte, name = expression(italic(P)-value~threshold)) + theme(axis.text.x=element_blank(), axis.title.x=element_blank())  +   xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  
    	            }
    	            if(bar.col.is.pval){
    	                ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = -log10(p.out)), stat="identity") +     scale_fill_gradient(low= bar.col.is.pval.lowcol, high= bar.col.is.pval.highcol, name =bquote(atop(-log[10]~model,italic(P)-value),)) +  xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  	        }
    	            ggfig.plot <- ggfig.plot + geom_text(aes(x = factor(thresh), y = r2.out, label = paste(print.p)), vjust = -1.5, hjust = 0, angle = 45, cex = 2.8, parse=T) +
    	                                        scale_y_continuous(limits = c(0, max(output$r2.out)*1.25))  + 
    	                                        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                axis.line = element_line(colour = "black",size=0.5))
                    if(binary.target){
                        if(!calculate.abc.r2){
                            ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
                        }
       	                if(calculate.abc.r2){
                            ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
                        }
                    }
                    if(!binary.target){
                        ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ",R^2)))
                    }
    	            ggfig.plot
    	            ggsave(paste(figname,"_",base.phenotypes.names[basePhen],"_BARPLOT_", Sys.Date(), ".png", sep = ""))		
    	        }  
    	        if(!ggfig) {
    	            png(paste(figname,"_",base.phenotypes.names[basePhen],"_BARPLOT_", Sys.Date(), ".png", sep = ""))
    	            plot.fig  <- with(output, barplot(r2.out, 
            	        names = thresh,
            	        main = "", 
            	        col = "red",
            	        xlab = expression(italic(P)[T]),
            	        ylab = expression(R^2),
            	        ylim = c(0,  max(output$r2.out)*1.25)  ))
            	        text( parse(text=paste(
            	        output$print.p)), 
            	        x = plot.fig+0.1, 
            	        y =  output$r2.out+ (max(output$r2.out)*1.125-max(output$r2.out)), 
            	        srt = 45)
    	            dev.off()
    	        }
    	    }
    	    if(multiple.target.phenotypes){
    	        for(k in 1:length(target.phenotypes)){
    	            output <- read.table(paste(figname, base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k], "RAW_RESULTS_DATA.txt", sep = "_"), head=T)
                    if(quantiles & ggfig){
                        cat(" ################################# \n # \n #   Quantiles Plots \n # \n ################################# \n")
                        if(binary.target){
                            if(report.individual.scores & !report.best.score.only){
                                scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
                            }
                            if(report.individual.scores & report.best.score.only){
                                scores.internal <- read.table(paste(figname,base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
                            }
                            names(scores.internal) <- c("ID", "SCORE")
                            for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
                            for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
                            for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))
    
                            if(covary == F){
                                or.quantiles <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1])
                                ci.quantiles.u <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
                                ci.quantiles.l <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
                            }
                            if(covary == T){
                                for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")
                                or.quantiles <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1])
                                ci.quantiles.u <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2]))
                                ci.quantiles.l <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2]))
                            }
                            or.quantiles[1] <- 1 
                            ci.quantiles.u[1] <- 1
                            ci.quantiles.l[1] <- 1
                            quant.list <- seq(1, num.quantiles, 1)
                            quant.list <- quant.list[quant.list != quant.ref]     
                            quantiles.for.table <- c(quant.ref, quant.list)
                            quantiles.df <- data.frame(or.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
                            names(quantiles.df) <- c("OR", "CI.U", "CI.L", "DEC")
                            quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
                            quantiles.plot <- ggplot(quantiles.df) + 
                              geom_point(aes(x = DEC, y = OR), colour = "royalblue2", size=4) + 
                              theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                              axis.line = element_line(colour = "black",size=0.5)) + 
                              geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = OR, x = DEC), colour = "royalblue2", size = 0.9) + 
                              ylab("Odds Ratio for Score on Phenotype") + 
                              xlab("quantiles for Polygenic Score") + 
                              scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
                        }  
                        if(!binary.target){
                            if(report.individual.scores & !report.best.score.only){
                                scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
                            }
                            if(report.individual.scores & report.best.score.only){
                                scores.internal <- read.table(paste(figname, base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
                            }
                            names(scores.internal) <- c("ID", "SCORE")
                            for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
                            for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
                            for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))
    
                            if(covary == F){
                                coef.quantiles <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1]
                                ci.quantiles.u <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
                                ci.quantiles.l <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
                            }
                            if(covary == T){
                                for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")
                                coef.quantiles <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] 
                                ci.quantiles.u <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
                                ci.quantiles.l <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
                            }
                            coef.quantiles[1] <- 0 
                            ci.quantiles.u[1] <- 0
                            ci.quantiles.l[1] <- 0
                            quant.list <- seq(1, num.quantiles, 1)
                            quant.list <- quant.list[quant.list != quant.ref]     
                            quantiles.for.table <- c(quant.ref, quant.list)
                            quantiles.df <- data.frame(coef.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
                            names(quantiles.df) <- c("Coef", "CI.U", "CI.L", "DEC")
                            quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
                            quantiles.plot <- ggplot(quantiles.df) + 
                              geom_point(aes(x = DEC, y = Coef), colour = "royalblue2", size=4) + 
                              theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                              axis.line = element_line(colour = "black",size=0.5)) + 
                              geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = Coef, x = DEC), colour = "royalblue2", size = 0.9) + 
                              ylab("Change in Phenotype given score in quantiles") + 
                              xlab("quantiles for Polygenic Score") + 
                              scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
                        }  
                        quantiles.plot
                        ggsave(paste(figname,base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k], "QUANTILES_PLOT.png", sep = "_"))  
                    }	
    
        	        if(!fastscore){
        	            cat(" ################################# \n # \n #   High Density Plots \n # \n ################################# \n")
        	        }
    	            if(best.thresh.on.bar){
            	        #add max thresh
            	        barchart.levels <- c(barchart.levels.old, output$thresh[which.min(output$p.out)])
            	        barchart.levels <- barchart.levels[!duplicated(barchart.levels)]	
            	        barchart.levels <- sort(barchart.levels, decreasing = F)
    	            }
    	            if(!fastscore & !ggfig){
    	                if(!scatter.R2){
                	        png(paste(figname,"_",base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
                	          with(output, plot(	x = output$thresh,
                	            y = -log10(output$p.out),
                	            xlab = expression(italic(P)-value~threshold),
                	            ylab = bquote(-log[10]~model~italic(P)-value),
                	            pch = 19,
                	            type = "b"
                	          ))
    	                    with(output[output$thresh %in% barchart.levels.old,], points(thresh, -log10(p.out), col = "green", pch = 19, type= "b"))
    	                }
    	                if(scatter.R2){
            	            png(paste(figname,"_",base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
            	            with(output, plot(	x = output$thresh,
                	            y = r2.out,
                	            xlab = expression(italic(P)-value~threshold),
                	            ylab = bquote(Variance~Explained:~Pseudo~R^2),
                	            pch = 19,
                	            type = "b"
            	            ))
    	                    with(output[output$thresh %in% barchart.levels.old,], points(thresh, r2.out, col = "green", pch = 19, type= "b"))
    	                }	
    	                dev.off()
    	            }  
    	            if(!fastscore & ggfig){
    	                if(!scatter.R2){
            	            ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = -log10(p.out))) +
                	            geom_line(aes(x = thresh, y = -log10(p.out))) +
                                xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
                                ylab(bquote(PRS~model~fit:~italic(P)-value~(-log[10]))) +
                	            geom_line(aes(thresh,  -log10(p.out)), colour = "green",
                	             data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
                	            geom_point(aes(thresh,  -log10(p.out)), colour = "green",
                	             data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
                	            theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                          axis.line = element_line(colour = "black",size=0.5))
            	            ggfig.points
            	            ggsave(paste(figname,"_",base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
    	                }
            	        if(scatter.R2){
            	            ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = r2.out)) +
                	            geom_line(aes(x = thresh, y = r2.out)) +
                                xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
                	            geom_line(aes(thresh,  r2.out), colour = "green",
                	            data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
                	            geom_point(aes(thresh,  r2.out), colour = "green",
                	            data = output[with(output, thresh %in% barchart.levels.old) , ] ) + 
                	            theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black",size=0.5))
                            if(binary.target){
                           	    if(!calculate.abc.r2){
                                    ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
                                }
                           	    if(calculate.abc.r2){
                                    ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
                                }
                            }
                            if(!binary.target){
                                ggfig.points <- ggfig.points + ylab(paste("PRS model fit:  ",expression(paste(R^2))))
                            }
    	                    ggfig.points
    	                    ggsave(paste(figname,"_",base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
    	                }	 
    	            }
    	            options(scipen=0,digits=7)
    	            output <- output[,c(1:3)]
    	            names(output) <- c("thresh", "p.out", "r2.out")
        	        output <- output[output$thresh %in% barchart.levels,]
        	        output$print.p[round(output$p.out, digits = 3) != 0] <- round(output$p.out[round(output$p.out, digits = 3) != 0], digits = 3)
        	        output$print.p[round(output$p.out, digits = 3) == 0] <- format(output$p.out[round(output$p.out, digits = 3) == 0], digits=2)
        	        output$print.p <- sub("e", "*x*10^", output$print.p)
        	        cat(" ################################# \n # \n #   Barplots \n # Bars for inclusion can be changed using the barchart.levels option \n ################################# \n")
    	            if(ggfig){
    	                if(!bar.col.is.pval){
    	                    ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = factor(thresh)), stat="identity") +     scale_fill_brewer(palette=barpalatte, name = expression(italic(P)-value~threshold)) + theme(axis.text.x=element_blank(), axis.title.x=element_blank())  +   xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  
    	                }
    	                if(bar.col.is.pval){
    	                    ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = -log10(p.out)), stat="identity") +     scale_fill_gradient(low= bar.col.is.pval.lowcol, high= bar.col.is.pval.highcol, name =bquote(atop(-log[10]~model,italic(P)-value),)) +   xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  
    	                }
    	                ggfig.plot <- ggfig.plot + geom_text(aes(x = factor(thresh), y = r2.out, label = paste(print.p)), vjust = -1.5, hjust = 0, angle = 45, cex = 2.8, parse=T) +
    	                                            scale_y_continuous(limits = c(0, max(output$r2.out)*1.25))  + 
                                        	        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                    axis.line = element_line(colour = "black",size=0.5))
                        if(binary.target){
       	                    if(!calculate.abc.r2){
                                ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
                            }
       	                    if(calculate.abc.r2){
                                ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
                            }
                        }
                        if(!binary.target){
                            ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ",R^2)))
                        }
    	                ggfig.plot
    	                ggsave(paste(figname,"_",base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],"_BARPLOT_", Sys.Date(), ".png", sep = ""))		
    	            }  
    	            if(!ggfig) {
        	            png(paste(figname,"_",base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],"_BARPLOT_", Sys.Date(), ".png", sep = ""))
        	            plot.fig  <- with(output, barplot(r2.out, 
            	            names = thresh,
            	            main = "", 
            	            col = "red",
            	            xlab = expression(italic(P)[T]),
            	            ylab = expression(R^2),
            	            ylim = c(0,  max(output$r2.out)*1.25)  ))
            	            text( parse(text=paste(
            	            output$print.p)), 
            	            x = plot.fig+0.1, 
            	            y =  output$r2.out+ (max(output$r2.out)*1.125-max(output$r2.out)), 
            	            srt = 45)
        	            dev.off()
    	            }
    	        }
    	    }
	    }
        if(multiple.base.phenotypes & (basePhen == 1) & !heat.r2){
            out.multibase <- data.frame(top.thresh.pval)
        }
        if(multiple.base.phenotypes & (basePhen > 1) & !heat.r2){
            out.multibase <- data.frame(out.multibase, top.thresh.pval)
        }
        if(multiple.base.phenotypes & (basePhen == 1) & heat.r2){
            out.multibase <- data.frame(top.thresh.r2)
        }
        if(multiple.base.phenotypes & (basePhen > 1) & heat.r2){
            out.multibase <- data.frame(out.multibase, top.thresh.r2)
        }
        if(multiple.base.phenotypes){
    	    cat(" ################################# \n # \n #   Remove temp files \n # \n ################################# \n")
    	    system("rm PROFILES.0.5*.profile", ignore.stdout=T,ignore.stderr=T)	
    	    system("rm PROFILES.0.4*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.0.3*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.0.2*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.0.1*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.0.0*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.S1*.profile", ignore.stdout=T,ignore.stderr=T)	
    	    system("rm PROFILES.S2*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.S3*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.S4*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.S5*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.S6*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.S7*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.S8*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.S9*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.S*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.0.5*.profile", ignore.stdout=T,ignore.stderr=T)	
    	    system("rm *PROFILES.0.4*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.0.3*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.0.2*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.0.1*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.0.0*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.S1*.profile", ignore.stdout=T,ignore.stderr=T)	
    	    system("rm *PROFILES.S2*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.S3*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.S4*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.S5*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.S6*.profile", ignore.stdout=T,ignore.stderr=T) 
    	    system("rm *PROFILES.S7*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.S8*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.S9*.profile", ignore.stdout=T,ignore.stderr=T)
    	    system("rm *PROFILES.S*.profile", ignore.stdout=T,ignore.stderr=T) 
    	    system("rm flip*", ignore.stdout=T,ignore.stderr=T)
    	    system("rm clean*", ignore.stdout=T,ignore.stderr=T)
    	    system("rm profile_list", ignore.stdout=T,ignore.stderr=T)
    	    system("rm rangelist.txt", ignore.stdout=T,ignore.stderr=T)
    	    system("rm rangelist_ranges", ignore.stdout=T,ignore.stderr=T)
    	    system("rm rawfile.raw", ignore.stdout=T,ignore.stderr=T)
    	    system("rm reordered_base", ignore.stdout=T,ignore.stderr=T)
    	    system("rm target_no_mhc*", ignore.stdout=T,ignore.stderr=T)
    	    system("rm ANCESTRY_INFORMATIVE_DIMENSIONS.cluster*", ignore.stdout=T,ignore.stderr=T)
    	    system("rm reordered_base", ignore.stdout=T,ignore.stderr=T)
    	    system("rm ANCESTRY_INFORMATIVE_DIMENSIONS.log", ignore.stdout=T,ignore.stderr=T)
    	    system("rm head_disc", ignore.stdout=T,ignore.stderr=T)
    	    system("rm LE_SNPs", ignore.stdout=T,ignore.stderr=T)
    	    system("rm prune_target*", ignore.stdout=T,ignore.stderr=T)
    	    system("rm plink.log", ignore.stdout=T,ignore.stderr=T)
    	    system("rm HEADER", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.hh", ignore.stdout=T,ignore.stderr=T)
    	    system("rm PROFILES.log", ignore.stdout=T,ignore.stderr=T)
    	    system("rm synonymous_snps*", ignore.stdout=T,ignore.stderr=T)
    	    system("rm synonymous_snps", ignore.stdout=T,ignore.stderr=T) 
    	    system("rm non_synonymous_snps_only*", ignore.stdout=T,ignore.stderr=T)
    	    system("rm mhc.txt", ignore.stdout=T,ignore.stderr=T)
    	    system("rm Rplots.pdf", ignore.stdout=T,ignore.stderr=T)
    	    system("rm temp.raw", ignore.stdout=T,ignore.stderr=T)
    	    system("rm TARGET_SNPs", ignore.stdout=T,ignore.stderr=T)
    	    system("rm base_SNPS", ignore.stdout=T,ignore.stderr=T)  
        }
    }
    if(multiple.base.phenotypes & multiple.target.phenotypes){
	    angle.heat <- 45
	    if(length(base.phenotypes.names) > 20){
		    angle.heat <- 90
	    }
        names(out.multibase) <- base.phenotypes.names
        row.names(out.multibase) <- target.phenotypes
        write.table(out.multibase, paste(figname, "ALL_BEST_THRESHOLDS_BASE_AND_TARGET.txt",sep="_"), col.names=T, row.names=T, quote= F)  
        if(ggfig){
            if(!heat.r2){
                trans.multibase <- data.frame(Base.Phenotype = rep(colnames(out.multibase), each = nrow(out.multibase)), Target.Phenotype = row.names(out.multibase), PRS.P.value = unlist(-log10(out.multibase)))
            }
            if(heat.r2){
                trans.multibase <- data.frame(Base.Phenotype = rep(colnames(out.multibase), each = nrow(out.multibase)), Target.Phenotype = row.names(out.multibase), PRS.P.value = unlist(out.multibase))
            }
            tile.plot <- ggplot(trans.multibase, aes(x = Base.Phenotype,
                                            y = Target.Phenotype, 
                                            fill = PRS.P.value)) + 
                                            geom_tile() +
                                            theme(axis.text.x = element_text(angle=angle.heat, hjust=1,vjust=1,size=12),
                                                axis.text.y = element_text(size=10),
                                                legend.text=element_text(size=10),
                                                legend.title=element_text(size=10)) 
            if(!heat.r2){
    	        tile.plot <- tile.plot + 
    	                        scale_fill_gradient(bquote(atop(Best~-log[10]~model,italic(P)-value),), 
        	                        low= bar.col.is.pval.lowcol, 
        	                        high= bar.col.is.pval.highcol)
            }
            if(heat.r2){
        	    tile.plot <- tile.plot + 
        	                scale_fill_gradient(bquote(Best~R^2), low= bar.col.is.pval.lowcol, high= bar.col.is.pval.highcol)
            }

            tile.plot
            ggsave(paste(figname, "_HEATMAP.png", sep=""))
        }
        if(!ggfig){
      	    png(paste(figname, "_HEATMAP.png", sep=""))
      	    heatmap(t(data.matrix(-log10(out.multibase))), Rowv=NA, Colv=NA, margins=c(10,10), cexCol=0.8, cexRow=0.8)
      	    dev.off()
        }
    }
    if(!(no.regression) & report.snps.best.model){
        all.raw <- read.table("rawfile.raw", head=F)
        names(all.raw) <- c("SNP", "A1", "BETA")
        all.p <- read.table("rangelist_ranges", head=F)
        names(all.p) <- c("SNP", "P")
        SNP.data <- merge(x = all.raw, y = all.p, by = "SNP")
        write.table(SNP.data, "ALL_SNPs_used_for_PRS.txt", col.names = T, row.names = F, quote=F)
        write.table(SNP.data[SNP.data$P < output$thresh[which.min(output$p.out)] , ], "SNPs_In_Best_model.txt",
            col.names = T, row.names = F, quote=F)
    }  
    if(cleanup){
        cat(" ################################# \n # \n #   Cleanup \n # \n ################################# \n")
        system("rm PROFILES.0.5*.profile", ignore.stdout=T,ignore.stderr=T)	
        system("rm PROFILES.0.4*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.0.3*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.0.2*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.0.1*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.0.0*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.S1*.profile", ignore.stdout=T,ignore.stderr=T)	
        system("rm PROFILES.S2*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.S3*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.S4*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.S5*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.S6*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.S7*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.S8*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.S9*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.S*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.0.5*.profile", ignore.stdout=T,ignore.stderr=T)	
        system("rm *PROFILES.0.4*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.0.3*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.0.2*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.0.1*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.0.0*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.S1*.profile", ignore.stdout=T,ignore.stderr=T)	
        system("rm *PROFILES.S2*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.S3*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.S4*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.S5*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.S6*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.S7*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.S8*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.S9*.profile", ignore.stdout=T,ignore.stderr=T)
        system("rm *PROFILES.S*.profile", ignore.stdout=T,ignore.stderr=T) 
        system("rm flip*", ignore.stdout=T,ignore.stderr=T)
        system("rm clean*", ignore.stdout=T,ignore.stderr=T)
        system("rm profile_list", ignore.stdout=T,ignore.stderr=T)
        system("rm rangelist.txt", ignore.stdout=T,ignore.stderr=T)
        system("rm rangelist_ranges", ignore.stdout=T,ignore.stderr=T)
        system("rm rawfile.raw", ignore.stdout=T,ignore.stderr=T)
        system("rm reordered_base", ignore.stdout=T,ignore.stderr=T)
        system("rm target_no_mhc*", ignore.stdout=T,ignore.stderr=T)
        system("rm ANCESTRY_INFORMATIVE_DIMENSIONS.cluster*", ignore.stdout=T,ignore.stderr=T)
        system("rm reordered_base", ignore.stdout=T,ignore.stderr=T)
        system("rm ANCESTRY_INFORMATIVE_DIMENSIONS.log", ignore.stdout=T,ignore.stderr=T)
        system("rm head_disc", ignore.stdout=T,ignore.stderr=T)
        system("rm LE_SNPs", ignore.stdout=T,ignore.stderr=T)
        system("rm prune_target*", ignore.stdout=T,ignore.stderr=T)
        system("rm plink.log", ignore.stdout=T,ignore.stderr=T)
        system("rm HEADER", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.hh", ignore.stdout=T,ignore.stderr=T)
        system("rm PROFILES.log", ignore.stdout=T,ignore.stderr=T)
        system("rm synonymous_snps*", ignore.stdout=T,ignore.stderr=T)
        system("rm synonymous_snps", ignore.stdout=T,ignore.stderr=T) 
        system("rm non_synonymous_snps_only*", ignore.stdout=T,ignore.stderr=T)
        system("rm mhc.txt", ignore.stdout=T,ignore.stderr=T)
        system("rm Rplots.pdf", ignore.stdout=T,ignore.stderr=T)
        system("rm temp.raw", ignore.stdout=T,ignore.stderr=T)
        system("rm TARGET_SNPs", ignore.stdout=T,ignore.stderr=T)
        system("rm base_SNPS", ignore.stdout=T,ignore.stderr=T)  
        system("rm Complete_Allele_List.txt", ignore.stdout=T,ignore.stderr=T)
    }
}

if(sumsum){
    cat(" ################################# \n # \n #   Using summary statistic - summary statistic method of Johnson et al \n # \n ################################# \n")

    if(plink.silent){
        plink <- paste(plink, " --silent ")
    }

    supper <- supper-sinc
    slower <- slower+sinc


    library(gtx)
    if(ggfig){
        library(ggplot2)  
        library(plyr)
    }
    options(stringsAsFactors=F)

    targ.gwas <- read.table(target,head=T)

    if(is.na(size.targ)){
        cat("ERROR: Please supply sample size for Target GWAS \n Quitting \n")
        quit()
    }

    if(clump.snps){
	    if(is.na(clump.ref)){
		    cat("ERROR: No reference panel supplied for clumping \n Quitting \n")
		    quit()
	    }
        system(paste("awk '{print $2}' ", clump.ref, ".bim | sort -k1,1 > clump_panel_snps.txt", sep = ""))
        system(paste("head -n +1 ", base, " > head.base"))
        system(paste("head -n +1 ", target, " > head.targ"))
        base.names <- read.table("head.base",head=T)
        targ.names <- read.table("head.targ",head=T)
        system(paste("tail -n +2 ", base, " | awk '{print $", which(colnames(base.names) == "SNP"), "}' | sort -k1,1 > base_snps",sep="")) 
        system(paste("tail -n +2 ", target, " | awk '{print $", which(colnames(targ.names) == "SNP"), "}' | sort -k1,1 > targ_snps",sep="")) 
        system(paste("join -1 1 -2 1 clump_panel_snps.txt base_snps | join -1 1 -2 1 ``-'' targ_snps > clean_snps", sep = ""))


        if("BETA" %in% colnames(base.names)){
            system(paste("tail -n +2 ", base, " | awk '{print $", which(colnames(base.names) == "SNP"),",$", which(colnames(base.names) == "A1"),",$", which(colnames(base.names) == "A2"),",$", which(colnames(base.names) == "BETA"),",$", which(colnames(base.names) == "P")," }' > reordered_base", sep=""))
            system("echo SNP A1 A2 BETA P > HEADER")
        }
        if(!("BETA" %in% colnames(base.names))){
            system(paste("tail -n +2 ", base, " | awk '{print $", which(colnames(base.names) == "SNP"),",$", which(colnames(base.names) == "A1"),",$", which(colnames(base.names) == "A2"),",$", which(colnames(base.names) == "OR"),",$", which(colnames(base.names) == "P")," }' > reordered_base", sep=""))
            system("echo SNP A1 A2 OR P > HEADER")
        }
        system(paste("sort -k1,1 reordered_base | join -1 1 -2 1 ``-'' clean_snps | cat HEADER ``-'' > cleaned_base", sep = ""))


        system(paste(plink," --noweb --bfile ", clump.ref, "  --clump cleaned_base --clump-p1 ",clump.p1," --clump-p2 ",clump.p2," --clump-r2 ",clump.r2," --clump-kb ",clump.kb," --out clumped_base", sep=" "))
        system("tail -n +2 clumped_base.clumped | awk '{print $3}'  | sort -k1,1 | awk '($1 != \"\"){print}'  > LE_SNPs")
        le.snps <- read.table("LE_SNPs",head=F)
        targ.gwas <- targ.gwas[targ.gwas$SNP %in% le.snps$V1 , ]
    }
    if(!clump.snps){
        system(paste("cp ", base, " cleaned_base"))
    }
    base.gwas <- read.table("cleaned_base",head=T)
    targ.gwas <- targ.gwas[targ.gwas$SNP %in% base.gwas$SNP , ]
    base.gwas <- base.gwas[base.gwas$SNP %in% targ.gwas$SNP , ]
    
    if("OR" %in% names(targ.gwas)){
        targ.gwas$BETA <- log(targ.gwas$OR)
    }
    if("OR" %in% names(base.gwas)){
        base.gwas$BETA <- log(base.gwas$OR)
    }

    cat(" ################################# \n # \n #   Check input format \n # \n ################################# \n")

    if(!("SNP" %in% names(base.gwas))){
        cat("ERROR: No SNP Name Column in Base GWAS \n Quitting \n")
        quit()
    }
    if(!("P" %in% names(base.gwas))){
        cat("ERROR: No P-value Column in Base GWAS \n Quitting \n")
        quit()
    }
    if(!("A1" %in% names(base.gwas))){
        cat("ERROR: No A1 Column in Base GWAS \n Quitting \n")
        quit()
    }
    if(!("A2" %in% names(base.gwas))){
        cat("ERROR: No A2 Column in Base GWAS \n Quitting \n")
        quit()
    }
    if(!("SNP" %in% names(targ.gwas))){
        cat("ERROR: No SNP Name Column in Target GWAS \n Quitting \n")
        quit()
    }
    if(!("SE" %in% names(targ.gwas))){
        cat("ERROR: No SE Column in Target GWAS \n Quitting \n")
    quit()
    }
    if(!("A1" %in% names(targ.gwas))){
        cat("ERROR: No A2 Column in Target GWAS \n Quitting \n")
        quit()
    }
    if(!("A2" %in% names(targ.gwas))){
        cat("ERROR: No A2 Column in Target GWAS \n Quitting \n")
        quit()
    }

    if(!("OR" %in% names(base.gwas) | "BETA" %in% names(base.gwas))){
        cat("ERROR: No Effect Size Column in Base GWAS \n ie OR or BETA \n Quitting \n")
        quit()
    }
    if(!("OR" %in% names(targ.gwas) | "BETA" %in% names(targ.gwas))){
        cat("ERROR: No Effect Size Column in Target GWAS \n ie OR or BETA \n Quitting \n")
        quit()
    }



    base.gwas <- base.gwas[order(base.gwas$SNP) , ]
    targ.gwas <- targ.gwas[order(targ.gwas$SNP) , ]

    base.gwas$BETA[base.gwas$A1 == targ.gwas$A2] <- -1*base.gwas$BETA[base.gwas$A1 == targ.gwas$A2]
    base.gwas <- base.gwas[base.gwas$A1 == targ.gwas$A1 | base.gwas$A1 == targ.gwas$A2  , ]


    high.res <- seq(slower, supper, sinc)
    pval.out <- as.vector(1)
    r2.out <- as.vector(1)
    nsnps <- as.vector(1)
    
    high.res <- c(high.res, barchart.levels)
    high.res <- sort(high.res)
    high.res <- high.res[!duplicated(high.res)]
    
    for(i in 1:length(high.res)){
        base.gwas.temp <- base.gwas[base.gwas$P < high.res[i] , ]
        targ.gwas.temp <- targ.gwas[targ.gwas$SNP %in% base.gwas.temp$SNP , ]
        pval.out[i] <- pnorm(abs(grs.summary(w = base.gwas.temp$BETA, b = targ.gwas.temp$BETA, s = targ.gwas.temp$SE, n = size.targ)$ahat/ grs.summary(w = base.gwas.temp$BETA, b = targ.gwas.temp$BETA, s = targ.gwas.temp$SE, n = size.targ)$aSE), lower.tail = F)
        r2.out[i] <- grs.summary(w = base.gwas.temp$BETA, b = targ.gwas.temp$BETA, s = targ.gwas.temp$SE, n = size.targ)$R2rs
        nsnps[i] <- dim(base.gwas.temp)[1]
    }

    output <- data.frame(high.res, pval.out, r2.out, nsnps)
    output$pval.out[output$pval.out == 0 ] <- 2.2e-16

    
    temp.output <- output
    names(temp.output) <- c("thresh", "pval", "r2", "nsnps")
    temp.output <- temp.output[!duplicated(temp.output$thresh),]
    temp.output <- temp.output[!is.na(temp.output$pval),]
    write.table(temp.output, paste(figname, "RAW_RESULTS_DATA.txt", sep = "_"), col.names=T, row.names=F,quote=F)
    
    cat(" ################################# \n # \n #   High Density Plots \n # \n ################################# \n")

    if(best.thresh.on.bar){
        #add max thresh
        barchart.levels <- c(barchart.levels.old, output$thresh[which.min(output$p.out)])
        barchart.levels <- barchart.levels[!duplicated(barchart.levels)]	
        barchart.levels <- sort(barchart.levels, decreasing = F)
    }
    
    if(ggfig){
        ggfig.points <- ggplot(data = output) + geom_point(aes(x = high.res, y = -log10(pval.out))) +
    	                    geom_line(aes(x = high.res, y = -log10(pval.out))) +
                            xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
                            ylab(bquote(Evidence~For~Shared~Genetic~Architecture:~italic(P)-value~(-log[10]))) +
    	                    geom_line(aes(high.res,  -log10(pval.out)), colour = "green",
    	                    data = output[with(output, high.res %in% barchart.levels) , ] ) + 
                            geom_point(aes(high.res,  -log10(pval.out)), colour = "green",
    	                    data = output[with(output, high.res %in% barchart.levels) , ] ) + 
    	                    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                            axis.line = element_line(colour = "black",size=0.5))
    	                    ggfig.points
        ggsave(paste(figname,"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
    }
    if(!ggfig){
        png(paste(figname,"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	          with(temp.output, plot(	x = temp.output$thresh,
	          y = -log10(temp.output$pval),
	          xlab = expression(italic(P)-value~threshold),
              ylab = bquote(Evidence~For~Shared~Genetic~Architecture:~italic(P)-value~(-log[10])),
	          pch = 19,
	          type = "b"
	        ))
	    with(temp.output[temp.output$thresh %in% barchart.levels,], points(thresh, -log10(pval), col = "green", pch = 19, type= "b"))
    }

    cat(" ################################# \n # \n #   Barplots \n # Bars for inclusion can be changed using the barchart.levels option \n ################################# \n")

    barchart.levels <- c(barchart.levels.old, output$high.res[which.min(output$pval.out)])
    barchart.levels <- sort(barchart.levels)
    
    output <- output[!duplicated(output$high.res),]
    output <- output[!is.na(output$pval.out),]
    
    
    options(scipen=0,digits=7)
    output <- output[output$high.res %in% barchart.levels,]
    output$print.p[round(output$pval.out, digits = 3) != 0] <- round(output$pval.out[round(output$pval.out, digits = 3) != 0], digits = 3)
    output$print.p[round(output$pval.out, digits = 3) == 0] <- format(output$pval.out[round(output$pval.out, digits = 3) == 0], digits=2)
    output$print.p <- sub("e", "*x*10^", output$print.p)
    if(ggfig){
        ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(high.res), y = r2.out, fill = -log10(pval.out)), stat="identity") +     
                                              scale_fill_gradient(low= bar.col.is.pval.lowcol, high= bar.col.is.pval.highcol, name =bquote(atop(-log[10]~model,italic(P)-value),)) +  
                                              xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  
                                              ggfig.plot <- ggfig.plot + geom_text(aes(x = factor(high.res), y = r2.out, label = paste(print.p)), vjust = -1.5, hjust = 0, angle = 45, cex = 2.8, parse=T) +
                                              scale_y_continuous(limits = c(0, max(output$r2.out)*1.25))  + 
                                            	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                      axis.line = element_line(colour = "black",size=0.5))
        if(binary.target){
	        ggfig.plot <- ggfig.plot + ylab(bquote(Variance~Explained:~Pseudo~R^2)) 
        } 
        if(binary.target){
	        ggfig.plot <- ggfig.plot + ylab(bquote(Variance~Explained:~R^2)) 
        } 
        ggfig.plot
        ggsave(paste(figname,"_BARPLOT_", Sys.Date(), ".png", sep = ""))		
    }
    if(!ggfig){
        png(paste(figname,"_BARPLOT_", Sys.Date(), ".png", sep = ""))
        if(binary.target){
            plot.fig  <- with(output, barplot(r2.out, 
                names = high.res,
                main = "", 
                col = "red",
                xlab = expression(italic(P)[T]),
                ylab = bquote(Pseudo~R^2),
                ylim = c(0,  max(output$r2.out)*1.25)  ))
                text( parse(text=paste(
                output$print.p)), 
                x = plot.fig+0.1, 
                y =  output$r2.out+ (max(output$r2.out)*1.125-max(output$r2.out)), 
                srt = 45)
        }
        if(!binary.target){
            plot.fig  <- with(output, barplot(r2.out, 
                            names = high.res,
                            main = "", 
                            col = "red",
                            xlab = expression(italic(P)[T]),
                            ylab = expression(R^2),
                            ylim = c(0,  max(output$r2.out)*1.25)  ))
                            text( parse(text=paste(
                            output$print.p)), 
                            x = plot.fig+0.1, 
                            y =  output$r2.out+ (max(output$r2.out)*1.125-max(output$r2.out)), 
                            srt = 45)
        }
        dev.off()
    }
    
    if(cleanup){
    	system("rm cleaned_base_temp", ignore.stdout=T,ignore.stderr=T)
    	system("rm clump_panel_snps.txt", ignore.stdout=T,ignore.stderr=T)
    	system("rm clumped_base.*", ignore.stdout=T,ignore.stderr=T)
    	system("rm head.base", ignore.stdout=T,ignore.stderr=T)
    	system("rm head.targ", ignore.stdout=T,ignore.stderr=T)
    	system("rm targ_snps", ignore.stdout=T,ignore.stderr=T)
    	system("rm HEADER", ignore.stdout=T,ignore.stderr=T)
    	system("rm cleaned_base", ignore.stdout=T,ignore.stderr=T)
    	system("rm clean_snps", ignore.stdout=T,ignore.stderr=T)
    	system("rm base_snps", ignore.stdout=T,ignore.stderr=T)
    	system("rm LE_SNPs", ignore.stdout=T,ignore.stderr=T)
    	system("rm Complete_Allele_List.txt", ignore.stdout=T,ignore.stderr=T)
    }
    	

}


if(print.time){
  cat(" ################################# \n # \n #   Print time \n # \n ################################# \n")
}

running.time <- proc.time()[3] - start.time 
if(running.time < 60){
	out.run.time <- round(running.time, digits = 2)
	out.time.units <- "seconds"
}
if(running.time > 60 & running.time < 3600){
	out.run.time <- round((running.time/60), digits = 2)
	out.time.units <- "minutes"
}
if(running.time > 3600){
	out.run.time <- round((running.time/3600), digits = 2)
	out.time.units <- "hours"
}


if(print.time){
  cat(paste("RUNNING TIME: ", 	out.run.time, out.time.units, "\n", sep = " "))
}



