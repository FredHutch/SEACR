#!/usr/bin/env bash

set -ue

display_help()
{
	echo "
	SEACR: Sparse Enrichment Analysis for CUT&RUN
	
	Usage: bash SEACR_1.4.sh -b <experimental bedgraph>.bg -c [<control bedgraph>.bg | <FDR threshold>] -n ["norm" | "non"] -m ["relaxed" | "stringent"] -o output_prefix -e [double] -r ["yes" | "no"]
	
		Required input fields:

		-b|--bedgraph Experimental bedgraph file
		-c|--control Control bedgraph file
		-n|--normalize Internal normalization (norm|non)
		-m|--mode Stringency mode (stringent|relaxed)
		-o|--output Output prefix

		Optional input fields:

		-e|--extentsion Peak extension constant for peak merging (Default=0.1)
		-r|--remove Remove peaks overlapping IgG peaks (yes|no, Default=yes)
		-h|--help	Print help screen
	
	Output file:
	<output prefix>.[stringent | relaxed].bed (Bed file of enriched regions)
	
	Output data structure: 
	
	<chr>	<start>	<end>	<AUC>	<max signal>	<max signal region>
	
	Description of output fields:
	Field 1: Chromosome
	
	Field 2: Start coordinate
	
	Field 3: End coordinate
	
	Field 4: Total signal contained within denoted coordinates
	
	Field 5: Maximum bedgraph signal attained at any base pair within denoted coordinates
	
	Field 6: Region representing the farthest upstream and farthest downstream bases within the denoted coordinates that are represented by the maximum bedgraph signal
	
	Examples:
	bash SEACR_1.4.sh -b target.bedgraph -c IgG.bedgraph -n norm -m stringent -o output
	Calls enriched regions in target data using normalized IgG control track with stringent threshold and default 10% peak extension for merging
	
	bash SEACR_1.4.sh -b target.bedgraph -c IgG.bedgraph -n non -m relaxed -o output -e 0.5
	Calls enriched regions in target data using non-normalized IgG control track with relaxed threshold and 50% peak extension for merging
	bash SEACR_1.4.sh -b target.bedgraph -c 0.01 -n non -m stringent -o output
	Calls enriched regions in target data by selecting the top 1% of regions by total signal
	"
	exit 1
}


POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -b|--bedgraph)
    BEDGRAPH="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--control)
   	CONTROL="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--normalize)
    NORMALIZE="$2"
    shift # past argument
    shift # past value
    ;;
    -m|--mode)
    MODE="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -e|--extension)
    EXTENSION="$2"
    shift # past argument
    shift # past value
    ;;
    -r|--remove)
    REMOVE="$2"
    shift # past argument
    shift # past value
    ;;
    -h|--help)
		display_help
		exit 1
    ;;
    --default)
    DEFAULT=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [[ -z ${EXTENSION+x} ]] 
then
	EXTENSION=0.1
fi

if [[ -z ${BEDGRAPH+x} ]] || [[ -z ${CONTROL+x} ]] || [[ -z ${NORMALIZE+x} ]] || [[ -z ${MODE+x} ]] || [[ -z ${OUTPUT+x} ]]
then
	echo "
		Missing required inputs!

		Required input fields:

		-b|--bedgraph Experimental bedgraph file
		-c|--control Control bedgraph file
		-n|--normalize Internal normalization (norm|non)
		-m|--mode Stringency mode (stringent|relaxed)
		-o|--output Output prefix

		Optional input fields:

		-e|--extentsion Peak extension constant for peak merging (Default=0.1)
		-r|--remove Remove peaks overlapping IgG peaks (yes|no, Default=yes)	
		-h|--help	Print help screen

	"
	exit 1
fi

if [[ -z ${REMOVE+x} ]]
then 
	REMOVE="yes"
fi
	
if [[ $BEDGRAPH == "-" ]]
then
	BEDGRAPH=$(cat)
fi

if [[ $REMOVE == "yes" ]] || [[ $REMOVE == "no" ]]
then
	echo "Removing peaks that overlap IgG peaks? $REMOVE"
else
	echo "--remove entry must be "yes" or "no"!"
	exit 1
fi

if [[ $EXTENSION =~ ^[0-9]?+([.][0-9]+)?$ ]] 
then
	echo "Extending peaks by $EXTENSION for merging"
else
	echo "--extension entry must be a positive number!"
	exit 1
fi

password=`head /dev/urandom | LC_CTYPE=C tr -dc A-Za-z0-9 | head -c 13; echo ''`
password2=`head /dev/urandom | LC_CTYPE=C tr -dc A-Za-z0-9 | head -c 13; echo ''`

exp=`basename $BEDGRAPH`

if [[ $CONTROL =~ ^[0-9]?+([.][0-9]+)?$ ]] || [[ $CONTROL =~ ^[0-9]([.][0-9]+) ]] || [[ $CONTROL =~ ^([.][0-9]+) ]]
then
	echo "Calling enriched regions without control file"
elif [[ -f $CONTROL ]]
then
	echo "Calling enriched regions with control file"
	ctrl=`basename $CONTROL`
else
	echo "$2 is not a number or a file"
	exit 1
fi

norm=`echo $NORMALIZE`

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

height=`echo $MODE`

if [[ $height == "relaxed" ]]
then
	echo "Using relaxed threshold"
elif [[ $height == "stringent" ]]
	then
	echo "Using stringent threshold"
else
	echo "Must specify \"stringent\" or \"relaxed\" in fourth input"
	exit 1
fi

echo "Creating experimental AUC file: $(date)"

awk 'BEGIN{s=1}; {if(s==1){s++}else if(s==2){if($4 > 0){chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2); num=1; s++}}else{if($4 > 0){if(chr==$1 && $2==stop){num++; stop=$3; auc=auc+($4*($3-$2)); if ($4 > max){max=$4; coord=$1":"$2"-"$3}else if($4 == max){split(coord,t,"-"); coord=t[1]"-"$3}}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord"\t"num; chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2); num=1}}}}' $BEDGRAPH > $password.auc.bed
cut -f 4,7 $password.auc.bed > $password.auc

if [[ -f $CONTROL ]]
then
  echo "Creating control AUC file: $(date)"

  awk 'BEGIN{s=1}; {if(s==1){s++}else if(s==2){if($4 > 0){chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2); num=1; s++}}else{if($4 > 0){if(chr==$1 && $2==stop){num++; stop=$3; auc=auc+($4*($3-$2)); if ($4 > max){max=$4; coord=$1":"$2"-"$3}else if($4 == max){split(coord,t,"-"); coord=t[1]"-"$3}}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord"\t"num; chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2); num=1}}}}' $CONTROL > $password2.auc.bed
  cut -f 4,7 $password2.auc.bed > $password2.auc
fi

# module load R  ## For use on cluster

echo "Calculating optimal AUC threshold: $(date)"

path=`dirname $0`
if [[ -f $CONTROL ]] && [[ $norm == "norm" ]]
then
	echo "Calculating threshold using normalized control: $(date)"
	Rscript $path/SEACR_1.4.R --exp=$password.auc --ctrl=$password2.auc --norm=yes --output=$password
elif [[ -f $CONTROL ]]
then
	echo "Calculating threshold using non-normalized control: $(date)"
	Rscript $path/SEACR_1.4.R --exp=$password.auc --ctrl=$password2.auc --norm=no --output=$password
else
	echo "Using user-provided threshold: $(date)"
	Rscript $path/SEACR_1.4.R --exp=$password.auc --ctrl=$2 --norm=no --output=$password
fi
	
fdr=`cat $password.fdr.txt | sed -n '1p'`			## Added 5/15/19 for SEACR_1.1
fdr2=`cat $password.fdr.txt | sed -n '2p'`			## Added 5/15/19 for SEACR_1.1

#thresh=`cat $exp.threshold.txt`
thresh=`cat $password.threshold.txt | sed -n '1p'`
thresh2=`cat $password.threshold.txt | sed -n '2p'`
thresh3=`cat $password.threshold.txt | sed -n '3p'`

echo "Creating thresholded feature file: $(date)"

if [[ $height == "relaxed" ]]
then
  echo "Empirical false discovery rate = $fdr2"
  awk -v value=$thresh2 -v value2=$thresh3 '$4 > value && $7 > value2 {print $0}' $password.auc.bed | cut -f 1,2,3,4,5,6 > $password.auc.threshold.bed
else
  echo "Empirical false discovery rate = $fdr"
  awk -v value=$thresh -v value2=$thresh3 '$4 > value && $7 > value2 {print $0}' $password.auc.bed | cut -f 1,2,3,4,5,6 > $password.auc.threshold.bed
fi

if [[ -f $CONTROL ]]
then
	if [[ $norm == "norm" ]] #If normalizing, multiply control bedgraph by normalization constant
	then
		constant=`cat $password.norm.txt | sed -n '1p'`
		awk -v mult=$constant 'BEGIN{OFS="\t"}; {$4=$4*mult; print $0}' $password2.auc.bed | cut -f 1,2,3,4,5,6 > $password2.auc2.bed
		mv $password2.auc2.bed $password2.auc.bed
	fi
	awk -v value=$thresh '$4 > value {print $0}' $password2.auc.bed > $password2.auc.threshold.bed
fi

echo "Merging nearby features and eliminating control-enriched features: $(date)"

# module load bedtools ## For use on cluster
linecount=`wc -l $password.auc.threshold.bed | awk '{print $1}'`
if [[ $linecount -gt 0 ]]
then
	mean=`awk -v scale=$EXTENSION '{s+=$3-$2; t++}END{print (s*scale)/t}' $password.auc.threshold.bed`
else
	mean=0
fi

if [[ -f $CONTROL ]] && [[ $REMOVE == "yes" ]]
then
	awk -v value=$mean 'BEGIN{s=1}; {if(s==1){chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6; s++}else{if(chr==$1 && $2 < stop+value){stop=$3; auc=auc+$4; if($5 > max){max=$5; coord=$6}else if($5==max){split(coord,t,"-"); split($6,u,"-"); coord=t[1]"-"u[2]}}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord; chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6}}}' $password.auc.threshold.bed | bedtools intersect -wa -v -a - -b $password2.auc.threshold.bed > $OUTPUT.auc.threshold.merge.bed  
else
	awk -v value=$mean 'BEGIN{s=1}; {if(s==1){chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6; s++}else{if(chr==$1 && $2 < stop+value){stop=$3; auc=auc+$4; if($5 > max){max=$5; coord=$6}else if($5==max){split(coord,t,"-"); split($6,u,"-"); coord=t[1]"-"u[2]}}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord; chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6}}}' $password.auc.threshold.bed > $OUTPUT.auc.threshold.merge.bed
fi

if [[ $height == "relaxed" ]]
then
  cat $OUTPUT.auc.threshold.merge.bed > $OUTPUT.relaxed.bed
else
  cat $OUTPUT.auc.threshold.merge.bed > $OUTPUT.stringent.bed
fi

echo "Removing temporary files: $(date)"

rm $password.auc.bed
rm $password.auc
rm $password.threshold.txt
rm $password.auc.threshold.bed
rm $password.fdr.txt  ## Added 5/15/19 for SEACR_1.1
rm $OUTPUT.auc.threshold.merge.bed
if [[ -f $CONTROL ]]
then
	rm $password2.auc.bed
	rm $password2.auc
	rm $password2.auc.threshold.bed
fi
if [[ $norm == "norm" ]]
then
	rm -f $password.norm.txt
fi
echo "Done: $(date)"
