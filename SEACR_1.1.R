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
invis <- gc(verbose=FALSE) 

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
suppressWarnings(numtest<-as.numeric(argsL$ctrl))
invis <- gc(verbose=FALSE)
if(is.na(numtest)){ ## If 2nd field is a bedgraph, calculate empirical threshold
#	print("Ctrl is a file")
	ctrl<-read.table(argsL$ctrl)
	ctrlvec<-ctrl$V1
	if(argsL$norm=="yes"){  ## Calculate peaks of density plots to generate normalization factor
		y<-seq(0,max(ctrlvec),length.out=max(ctrlvec)-1)
		ctrlvalue<-max(which(ecdf(ctrlvec)(y)<=0.9))
		expvalue<-max(which(ecdf(expvec)(y)<=0.9))
		ctrltest<-density(ctrlvec[ctrlvec <= ctrlvalue]) ## New for SEACR_1.1
		exptest<-density(expvec[expvec <= expvalue]) ## New for SEACR_1.1
		constant<-(exptest$x[exptest$y==max(exptest$y)])/(ctrltest$x[ctrltest$y==max(ctrltest$y)])
		ctrlvec<-ctrlvec*constant
	} ## Calculate total signal and max signal thresholds
	both<-c(expvec,ctrlvec)
	pctremain<-function(x) (length(expvec)-(ecdf(expvec)(x)*length(expvec)))/(length(both)-(ecdf(both)(x)*length(both)))
	x<-sort(unique(both)) ## New for SEACR_1.1
	x0<-x[which(na.omit(pctremain(x[pctremain(x) < 1])) == max(na.omit(pctremain(x[pctremain(x) < 1]))))]  ## New for SEACR_1.1
	z<-x[x <= x0[1]]	## New for SEACR_1.1
	z2<-z[abs(((pctremain(x0)+min(pctremain(z)))/2)-pctremain(z))==min(abs(((pctremain(x0)+min(pctremain(z)))/2)-pctremain(z)))] ## New for SEACR_1.1
	z<-z[z > z2[1]] ## New for SEACR_1.1
	z0<-z[abs(z-(max(z)-((1/2)*(max(z)-min(z)))))==min(abs(z-(max(z)-((1/2)*(max(z)-min(z))))))] ## New for SEACR_1.1
	
	## The following code segment was added to avoid spurious high thresholding when the peak of a lower threshold is within 95% of the peak of the maximum threshold
	
	frame<-data.frame(thresh=x[1:(length(x)-1)], pct=pctremain(x[1:(length(x)-1)]), diff=abs(diff(pctremain(x))))
	frame<-na.omit(frame)
	i<-2
	output<-0
	while(output==0){
		test<-as.numeric(paste(c(0,".",rep(9,i)),sep="",collapse=""))
		output<-as.numeric(quantile(frame$diff, test))
#		print(output)
		i<-i+1
	}
	a<-frame$thresh[frame$diff != 0 & frame$diff < quantile(frame$diff, test)]
	a0<-a[which(na.omit(pctremain(a[pctremain(a) < 1])) == max(na.omit(pctremain(a[pctremain(a) <  1]))))]
	b<-a[a <= a0[1]]
	b2<-b[abs(((pctremain(a0)+min(pctremain(b)))/2)-pctremain(b))==min(abs(((pctremain(a0)+min(pctremain(b)))/2)-pctremain(b)))]
	b<-b[b > b2[1]]
	b0<-b[abs(b-(max(b)-((1/2)*(max(b)-min(b)))))==min(abs(b-(max(b)-((1/2)*(max(b)-min(b))))))]
	if(max(na.omit(pctremain(a[pctremain(a) < 1])))/max(na.omit(pctremain(x[pctremain(x) < 1]))) > 0.95){
		x0<-a0
		z0<-b0
	}
	invis <- gc(verbose=FALSE)
	fdr<-c(1-pctremain(x0[1]), 1-pctremain(z0[1])) ## New for SEACR_1.1
}else{ ## If 2nd field is numeric, calculate percentile threshold
#	print("Ctrl is numeric")
	test<-ecdf(exp$V1)(exp$V1)
	frame<-data.frame(values=exp$V1, percentile=1-test)
	test2<-ecdf(exp$V2)(exp$V2)
	frame2<-data.frame(values=exp$V2, percentile=1-test2)
	ctrl<-as.vector(argsL$ctrl)
	x0<-min(frame$values[frame$percentile <= ctrl[1]])
	z0<-min(frame2$values[frame2$percentile <= ctrl[1]])
	fdr<-ctrl[1] ## New for SEACR_1.1
}
invis <- gc(verbose=FALSE)
write.table(c(x0[1],z0[1]), file=paste(argsL$output, ".threshold.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
if(argsL$norm=="yes"){
	write.table(constant, file=paste(argsL$output, ".norm.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) #Added 7/19/18 to ensure norm value is multiplied by ctrl
}
invis <- gc(verbose=FALSE)
write.table(fdr, file=paste(argsL$output, ".fdr.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) #Added 5/15/19 to report empirical FDR for threshold detection
