#!/bin/bash
SAMPLE=$1
BAMFILE=$2
SOURCEBAM=$3
PROJECT=$4
LUMPYPATH=$5
BAMDIR=`dirname $BAMFILE`
LB=`sambamba view -H $BAMDIR/$SAMPLE.bam | grep -m 1 "^@RG" | awk '{ for (i=1;i<=NF;++i) { if ($i~"^LB:") LB=$i; gsub("^LB:","",LB); } print LB }'`
SM_PRINT=`sambamba view -H $BAMDIR/$SAMPLE.bam | grep -m 1 "^@RG" | awk '{ for (i=1;i<=NF;++i) { if ($i~"^SM:") SM=$i; gsub("^SM:","",SM); } print SM }'`
echo $LB $SM_PRINT
LUMPYPATH=$LUMPYPATH$SAMPLE
# 	#TRANSPOSE
flagstat=`awk '
{
for (i=1; i<=NF; i++)  {
a[NR,i] = $i
}
}
NF>p { p = NF }
END {
for(j=1; j<=p; j++) {
str=a[1,j]
for(i=2; i<=NR; i++){
str=str" "a[i,j];
}
print str
}
}' $BAMDIR/$SAMPLE.flagstat | head -n 1 | awk -v var=$SAMPLE 'BEGIN {OFS="\t"} {print $0 }'`
	
	#total	duplicates	mapped	paired	read1	read2	properly-paired	mate-mapped	singletons	mate.diff.chr	mate.diff.chr.mapQ5
	echo $flagstat
	# #REST OF CALCULATIONS
	genomelen=3101976562
	export primarymapped=`echo $flagstat | cut -d ' ' -f 15`
	export meaninssz=`echo $flagstat | cut -d ' ' -f 16`
	export sdinssz=`echo $flagstat | cut -d ' ' -f 17`
	export meanedit=`echo $flagstat | cut -d ' ' -f 18`
	export readlen=`cat $BAMDIR/$SAMPLE.readlen | cut -f 1 | tail -n 1 | awk '{if (NR == 1) print $1;}'`
	genomecov=`R -q -e  "cat($primarymapped/$genomelen*$readlen)"`
	genomecov=`echo $genomecov | sed 's/cat([^)]*)\s//g;s/>//g'`
	physcov=`R -q -e "cat($primarymapped/2*$meaninssz/$genomelen)"`
	physcov=`echo $physcov | sed 's/cat([^)]*)\s//g;s/>//g'`
	errorrate=`R -q -e "cat($meanedit/$readlen)"`
	errorrate=`echo $errorrate | sed 's/cat([^)]*)\s//g;s/>//g'`
	#echo -e $SM_PRINT"\t"$LB"\t"$BAMDIR"\t"$SOURCEBAM"\t"$PROJECT"\t"$LUMPYPATH"\t"$flagstat"\t"$readlen"\t"$genomecov"\t"$physcov"\t"$errorrate  >> stats.summary
	echo -e $SM_PRINT"\t"$LB"\t"$BAMDIR"/"$SAMPLE"/"$SAMPLE".bam\t"$SOURCEBAM"\t"$PROJECT"\t"$LUMPYPATH"\t"$flagstat"\t"$readlen"\t"$genomecov"\t"$physcov"\t"$errorrate | tr ' ' '\t' | awk '{$8=$9=$20="";print $0}' >>stats.new
	