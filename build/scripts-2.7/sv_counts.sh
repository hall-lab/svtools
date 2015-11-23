#!/bin/bash

VCF=$1
SAMPLE=$2
QUAL=$3
SM_PRINT=`basename $VCF`
SM_PRINT=${SM_PRINT%.*.*.*.*}
if [[ -z "$VCF" ]] || [[ -z "$SAMPLE" ]]
then
    echo "usage: $0 <VCF> <SAMPLE>"
    exit 1
fi
if [[${VCF##*.} == ".gz"]]
then
	for TYPE in DEL DUP INV BND
	do
	zcat $VCF \
	| vawk -v TYPE=$TYPE -v QUAL=$QUAL "{ if (I\$SVTYPE==TYPE && ! I\$SECONDARY && \$6>=QUAL) print S\$$SAMPLE\$GT}"  \
	| sort | uniq -c \
	| awk -v SLOP=$SLOP -v SAMPLE=$SM_PRINT -v TYPE=$TYPE 'BEGIN {HOMREF=0; HET=0; HOMALT=0} { if ($2=="0/0") HOMREF=$1; else if ($2=="0/1") HET=$1; else if ($2=="1/1") HOMALT=$1 } END {print SAMPLE,SLOP,TYPE,HOMREF,HET,HOMALT}' OFS="\t"
	done
else
	for TYPE in DEL DUP INV BND
	do
	cat $VCF \
	| vawk -v TYPE=$TYPE -v QUAL=$QUAL "{ if (I\$SVTYPE==TYPE && ! I\$SECONDARY && \$6>=QUAL) print S\$$SAMPLE\$GT}"  \
	| sort | uniq -c \
	| awk -v SLOP=$SLOP -v SAMPLE=$SM_PRINT -v TYPE=$TYPE 'BEGIN {HOMREF=0; HET=0; HOMALT=0} { if ($2=="0/0") HOMREF=$1; else if ($2=="0/1") HET=$1; else if ($2=="1/1") HOMALT=$1 } END {print SAMPLE,SLOP,TYPE,HOMREF,HET,HOMALT}' OFS="\t"
	done
fi