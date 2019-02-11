#!/usr/bin/bash

set -ue

if [ $# -lt 4 ]
then
	echo "
	SEARCh: Sparse Enrichment Analysis for Regions in Chromatin
	
	Usage: bash SEARCh_1.0.sh <experimental bedgraph>.bg [<control bedgraph>.bg | <FDR threshold>] ["norm" | "non"] ["union" | "AUC"]
	
	Output:

	<experimental bedgraph>.auc.threshold.merge.bed (Bed file of enriched regions)
	Data structure: <chr>	<start>	<stop>	<AUC>	<max signal>	<max signal region>
	
	<experimental bedgraph>.auc.threshold.merge.summits.bed (Bed file of single base-pair summits)
	Data structure: <chr>	<start>	<stop>	<Predicted fragment length>	<Predicted dispersion (sd/mean)>

	Example:

	bash SEARCh_1.0.sh MPM091.bedgraph MPM099.bedgraph norm AUC
	Calls enriched regions in MPM091 using normalized IgG control track from MPM099
	
	bash SEARCh_1.0.sh MPM091.bedgraph MPM099.bedgraph non union
	Calls enriched regions in MPM091 using non-normalized IgG control track from MPM099 

	bash SEARCh_1.0.sh MPM091.bedgraph 0.01 non AUC
	Calls enriched regions in MPM091 by selecting the top 1% of regions by area under the curve (AUC)

	"
	exit 1
fi

password=`head /dev/urandom | tr -dc A-Za-z0-9 | head -c 13; echo ''`
password2=`head /dev/urandom | tr -dc A-Za-z0-9 | head -c 13; echo ''`

exp=`basename $1`

if [[ $2 =~ ^[0-9]+([.][0-9]+)?$ ]]
then
	echo "Calling enriched regions without control file"
elif [[ -f $2 ]]
then
	echo "Calling enriched regions with control file"
	ctrl=`basename $2`
else
	echo "$2 is not a number or a file"
	exit 1
fi

if [ $# -eq 5 ]
then
	echo "Using fragment size-based summit detection"
fi

norm=`echo $3`

if [[ $norm == "norm" ]]
then
	echo "Normalizing control to experimental bedgraph"
elif [[ $norm == "non" ]]
	then
	echo "Proceeding without normalization of control to experimental bedgraph"
else
	echo "Must specify \"norm\" for normalized or \"non\" for non-normalized data processing in third input"
	exit 1
fi

height=`echo $4`

if [[ $height == "union" ]]
then
	echo "Using peak height in addition to AUC threshold"
elif [[ $height == "AUC" ]]
	then
	echo "Proceeding without peak height threshold"
else
	echo "Must specify \"union\" for peak height threshold or \"AUC\" for no peak height threshold in fourth input"
	exit 1
fi

echo "Creating experimental AUC file: $(date)"

awk 'BEGIN{s=1}; {if(s==1){s++}else if(s==2){chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2); s++}else{if(chr==$1 && $2==stop){stop=$3; auc=auc+($4*($3-$2)); if ($4 > max){max=$4; coord=$1":"$2"-"$3}else if($4 == max){split(coord,t,"-"); coord=t[1]"-"$3}}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord; chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2)}}}' $1 > $password.auc.bed
cut -f 4,5 $password.auc.bed > $password.auc

if [[ -f $2 ]]
then
	echo "Creating control AUC file: $(date)"

	awk 'BEGIN{s=1}; {if(s==1){s++}else if(s==2){chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2); s++}else{if(chr==$1 && $2==stop){stop=$3; auc=auc+($4*($3-$2)); if ($4 > max){max=$4; coord=$1":"$2"-"$3}else if($4 == max){split(coord,t,"-"); coord=t[1]"-"$3}}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord; chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2)}}}' $2 > $password2.auc.bed
	cut -f 4,5 $password2.auc.bed > $password2.auc
fi

module load R

echo "Calculating optimal AUC threshold: $(date)"

if [[ -f $2 ]] && [[ $norm == "norm" ]]
then
	echo "Calculating threshold using normalized control: $(date)"
	Rscript SEARCh_1.0.R --exp=$password.auc --ctrl=$password2.auc --norm=yes --output=$password
elif [[ -f $2 ]]
then
	echo "Calculating threshold using non-normalized control: $(date)"
	Rscript SEARCh_1.0.R --exp=$password.auc --ctrl=$password2.auc --norm=no --output=$password
else
	echo "Using user-provided threshold: $(date)"
	Rscript SEARCh_1.0.R --exp=$password.auc --ctrl=$2 --norm=no --output=$password
fi
	
#thresh=`cat $exp.threshold.txt`
thresh=`cat $password.threshold.txt | sed -n '1p'`
thresh2=`cat $password.threshold.txt | sed -n '2p'`

echo "Creating thresholded feature file: $(date)"

if [[ $height == "union" ]] 
then
#	awk -v value=$thresh '$4 > value {print $0}' $password.auc.bed | awk -v value2=$thresh2 '$5 < value2 {print $0}' > $password.auc.threshold.bed (Previous behavior)
	awk -v value=$thresh -v value2=$thresh2 '$4 > value || $5 > value2 {print $0}' $password.auc.bed > $password.auc.threshold.bed #(Current behavior as of 2/7/19)
else
	awk -v value=$thresh '$4 > value {print $0}' $password.auc.bed > $password.auc.threshold.bed 
fi

if [[ -f $2 ]]
then
#	if [[ $height == "union" ]] 
#	then
#		awk -v value=$thresh -v value2=$thresh2 '$4 > value || $5 > value2 {print $0}' $password2.auc.bed > $password2.auc.threshold.bed
#	else
	if [[ $norm == "norm" ]] #If normalizing, multiply control bedgraph by normalization constant
	then
		constant=`cat $password.norm.txt | sed -n '1p'`
		awk -v mult=$constant 'BEGIN{OFS="\t"}; {$4=$4*mult; print $0}' $password2.auc.bed > $password2.auc2.bed
		mv $password2.auc2.bed $password2.auc.bed
	fi
	awk -v value=$thresh '$4 > value {print $0}' $password2.auc.bed > $password2.auc.threshold.bed
#	fi
fi

echo "Merging nearby features and eliminating control-enriched features: $(date)"

module load bedtools
mean=`awk '{s+=$3-$2; t++}END{print s/(t*10)}' $password.auc.threshold.bed`

if [[ -f $2 ]]
then
	awk -v value=$mean 'BEGIN{s=1}; {if(s==1){chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6; s++}else{if(chr==$1 && $2 < stop+value){stop=$3; auc=auc+$4; if($5 > max){max=$5; coord=$6}else if($5==max){split(coord,t,"-"); split($6,u,"-"); coord=t[1]"-"}u[2]}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord; chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6}}}' $password.auc.threshold.bed | bedtools intersect -wa -v -a - -b $password2.auc.threshold.bed > $exp.auc.threshold.merge.bed  
else
	awk -v value=$mean 'BEGIN{s=1}; {if(s==1){chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6; s++}else{if(chr==$1 && $2 < stop+value){stop=$3; auc=auc+$4; if($5 > max){max=$5; coord=$6}else if($5==max){split(coord,t,"-"); split($6,u,"-"); coord=t[1]"-"}u[2]}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord; chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6}}}' $password.auc.threshold.bed > $exp.auc.threshold.merge.bed
fi

echo "Removing temporary files: $(date)"

rm $password.auc.bed
rm $password.auc
rm $password.threshold.txt
rm $password.auc.threshold.bed
if [[ -f $2 ]]
then
	rm $password2.auc.bed
	rm $password2.auc
	rm $password2.auc.threshold.bed
fi
if [[ $norm == "norm" ]]
then
	rm $password.norm.txt
fi
echo "Done: $(date)"