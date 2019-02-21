#!usr/bin/Rscript

## Collect arguments
args <- commandArgs(TRUE)
 
## Default setting when no arguments passed
if(length(args) < 4) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
     Calculate area under the curve threshold for CUT&RUN peaks 
 
      Arguments:
			--exp=someValue   - Input AUC values from experiment CUT&RUN
			--ctrl=someValue   - Input AUC values from control CUT&RUN
			--norm=[yes|no]     - Whether to normalize control and experimental files
			--output=someValue   - Output prefix
")
 
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
 
## Arg1 default
#if(is.null(args[1])){
if(is.null(argsL$exp) | is.null(argsL$ctrl) | is.null(argsL$output) | is.null(argsL$norm)) {
  stop("Argument is missing!
     Calculate area under the curve threshold for CUT&RUN peaks 
 
      Arguments:
			--exp=someValue   - Input AUC values from experiment CUT&RUN
			--ctrl=someValue   - Input AUC values from control CUT&RUN
			--norm=[yes|no]     - Whether to normalize control and experimental files
			--output=someValue   - Output prefix
")
 
  q(save="no")
}
exp<-read.table(argsL$exp)
expvec<-exp$V1
expmax<-exp$V2
numtest<-as.numeric(argsL$ctrl)
if(is.na(numtest)){ ## If 2nd field is a bedgraph, calculate empirical threshold
#	print("Ctrl is a file")
	ctrl<-read.table(argsL$ctrl)
	ctrlvec<-ctrl$V1
	ctrlmax<-ctrl$V2
	if(argsL$norm=="yes"){  ## Calculate peaks of density plots to generate normalization factor
		y<-seq(0,max(ctrlvec),length.out=max(ctrlvec)-1)
		ctrlvalue<-max(which(ecdf(ctrlvec)(y)<=0.9))
		expvalue<-max(which(ecdf(expvec)(y)<=0.9))
		value<-max(c(expvalue,ctrlvalue))
		ctrltest<-density(ctrlvec[ctrlvec <= value])
		exptest<-density(expvec[expvec <= value])
		constant<-(exptest$x[exptest$y==max(exptest$y)])/(ctrltest$x[ctrltest$y==max(ctrltest$y)])
		ctrlvec<-ctrlvec*constant
		ctrlmax<-ctrlmax*constant
	} ## Calculate total signal and max signal thresholds
	both<-c(expvec,ctrlvec)
	both2<-c(expmax,ctrlmax)
	pctremain<-function(x) (length(expvec)-(ecdf(expvec)(x)*length(expvec)))/(length(both)-(ecdf(both)(x)*length(both)))
	x<-seq(0,floor(max(ctrlvec))-1,length.out=floor(max(ctrlvec))-1)
	x0<-x[which(na.omit(pctremain(x)) == max(na.omit(pctremain(x))))]
	pctremain2<-function(x) (length(expmax)-(ecdf(expmax)(x)*length(expmax)))/(length(both2)-(ecdf(both2)(x)*length(both2)))
	if(max(ctrlmax) < 1){
		z<-seq(0,max(ctrlmax),length.out=1000)
	}else{
		z<-seq(0,floor(max(ctrlmax))-1,length.out=floor(max(ctrlmax))-1)
	}
	z0<-z[which(na.omit(pctremain2(z)) == max(na.omit(pctremain2(z))))]
}else{ ## If 2nd field is numeric, calculate percentile threshold
#	print("Ctrl is numeric")
	test<-ecdf(exp$V1)(exp$V1)
	frame<-data.frame(values=exp$V1, percentile=1-test)
	test2<-ecdf(exp$V2)(exp$V2)
	frame2<-data.frame(values=exp$V2, percentile=1-test2)
	ctrl<-as.vector(argsL$ctrl)
	x0<-min(frame$values[frame$percentile <= ctrl[1]])
	z0<-min(frame2$values[frame2$percentile <= ctrl[1]])
}
write.table(c(x0[1],z0[1]), file=paste(argsL$output, ".threshold.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
if(argsL$norm=="yes"){
	write.table(constant, file=paste(argsL$output, ".norm.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) #Added 7/19/18 to ensure norm value is multiplied by ctrl
}
