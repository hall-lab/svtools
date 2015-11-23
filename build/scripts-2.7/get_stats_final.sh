#! /gsc/bin/bash
SAMPLE=$1
DIR=$2
OP=$3
echo -e $SAMPLE "\n" $DIR "\n" $SCHEDULER "\n"
#insert_size
bomb -q hall-lab -m 2 -g /abadve/stats -J $SAMPLE"-insert" -o $OP$SAMPLE/$SAMPLE.log -e $OP$SAMPLE/$SAMPLE.log \
 	  	"sambamba view -F 'paired and mate_is_reverse_strand and not (unmapped or reverse_strandor duplicate or secondary_alignment or mate_is_unmapped)' $DIR$SAMPLE/$SAMPLE.bam | awk '{if (\$7==\"=\" && \$9>0 && \$9<1000) print \$9;}' | head -n 10000000 | tail -n 5000000 > $OP$SAMPLE/$SAMPLE.insert.size"
bomb -q hall-lab -m 4 -g /abadve/stats  -J $SAMPLE"-NM"  -o $OP$SAMPLE/$SAMPLE.log -e $OP$SAMPLE/$SAMPLE.log \
	"sambamba view -f bam -F 'paired and mate_is_reverse_strand and not (unmapped or reverse_strand or duplicate or secondary_alignment or mate_is_unmapped)' $DIR$SAMPLE/$SAMPLE.bam | bedtools bamtobed -i stdin -tag NM | cut -f 5 | head -n 1000000 | tail -n 5000000 > $OP$SAMPLE/$SAMPLE.edit.dist"
bomb -q hall-lab -m 1 -g /abadve/stats -J $SAMPLE"-flag" -o $OP$SAMPLE/$SAMPLE.log -e $OP$$SAMPLE/$SAMPLE.log \
		"sambamba flagstat $DIR$SAMPLE/$SAMPLE.bam > $OP$SAMPLE/$SAMPLE.flagstat"
# bomb -q hall-lab -m 2 -g /abadve/stats -J $SAMPLE"-flag" -o $OP/$SAMPLE/$SAMPLE.log -e $OP/$$SAMPLE/$SAMPLE.log \
# 	"/gscuser/habel/src/samtools-develop/samtools flagstat_mod $DIR/$SAMPLE/$SAMPLE.bam >$OP/$SAMPLE/$SAMPLE.flagstat"
#PRIMARY MAPPING
bomb -q hall-lab -m 2 -g /abadve/stats  -J $SAMPLE"-prim"  -o $OP$SAMPLE/$SAMPLE.log -e $OP$SAMPLE/$SAMPLE.log \
	"sambamba view -F 'not (unmapped or secondary_alignment)' $DIR$SAMPLE/$SAMPLE.bam | wc -l > $OP$SAMPLE/$SAMPLE.primarymapped"
	# #readlength
bomb -q hall-lab -g /abadve/stats  -J $SAMPLE"-readlen" -o $OP/$SAMPLE/$SAMPLE.log -e $OP/$SAMPLE/$SAMPLE.log \
	"sambamba view -t 3 $DIR/$SAMPLE/$SAMPLE.bam | perl -lane '\$l = 0; \$F[5] =~s/(\d+)[MSINHP]/\$l+=\$1/eg; print \$l' | head -n 1 > $OP/$SAMPLE/$SAMPLE.readlen"

# paired
# proper_pair
# unmapped
# mate_is_unmapped
# reverse_strand
# mate_is_reverse_strand
# first_of_pair
# second_of_pair
# secondary_alignment
# failed_quality_control
# duplicate
# supplementary
# Bit Description
# 0x1 template having multiple segments in sequencing
# 0x2 each segment properly aligned according to the aligner
# 0x4 segment unmapped
# 0x8 next segment in the template unmapped
# 0x10 SEQ being reverse complemented
# 0x20 SEQ of the next segment in the template being reversed
# 0x40 the first segment in the template
# 0x80 the last segment in the template
# 0x100 secondary alignment
# 0x200 not passing quality controls
# 0x400 PCR or optical duplicate
# 0x800 supplementary alignment