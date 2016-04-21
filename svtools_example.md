#Example Workflow on Illumina Platinum Data

## Install svtools using pip package from pypi
First you will need to prepare your python environment.
You might want to use pyenv virtualenv.
You can learn about that here <https://github.com/yyuu/pyenv-virtualenv>
We require python 2.7.X
The creation of the pyenv and activation looks like this in our environment
<pre><code>pyenv virtualenv 2.7.9 svtools_install_instructions-2.7.9
pyenv activate svtools_install_instructions-2.7.9</pre></code>
Now you will need to satisfy the pysam dependency
<pre><code>pip install pysam>=0.8.1,<0.9.0</code><pre>
Then you should be able to install the svtools package from pypi
<pre><code>pip install svtools</pre></code>
You can spot check your svtools install by running
<pre><code>svtools --version</code></pre>
additional installtion strategies available in DEVELOPER.md

##Download bamfiles for NA12878 pedigrees from http://www.ebi.ac.uk/ena/data/view/ERP001960
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/NA12877_S1.bam
that is gonna take a WHILE, a single bam looked to take about 18 hours and the total number of bams suggested by this instruction is 
about 14 days worth of downloads from machines like our workstations....perhaps quite a bit shorter from blades

## Create a bampaths mapping file with original paths to all unaligned NA12878 pedigree bam files

## Alternative directory of realigned BAM files at hall-lab disk to use for analysis
   # /gscmnt/gc2802/halllab/ceph1463_realign_021815/
   # /gscmnt/gc2802/halllab/ceph1463_realign_021815/notes/sample.path.txt

# (1) make the directories
for SAMPLE in `cat $SAMPLEMAP | awk '{print $1}'` ; do mkdir -p $WORKDIR/$SAMPLE/log; \
	 mkdir -p $WORKDIR/$SAMPLE/qc; done

#speedseq realign allows alignment from one or more BAM files, rather than FASTQ inputs. It automatically read group information from the BAM header to mark duplicates by library.

# ----------------------------------------
# 1. speedseq
# ----------------------------------------
# Program: speedseq
# Local Path : /gscmnt/gc2719/halllab/bin/speedseq
# Version: 0.0.3
# Author: Colby Chiang
# usage:   speedseq <command> [options]
# command: align    align FASTQ files with BWA-MEM 
#          var      call SNV and indel variants with FreeBayes
#          somatic  call somatic SNV and indel variants in a tumor/normal pair with FreeBayes
#          sv       call SVs with LUMPY
#          realign  realign from a coordinate sorted BAM file 
#
# Dependencies on other software/packages:
# 1. BWA: speedseq uses bwa mem to   
	# Program   : bwa mem:
	# Local Path			: /gscmnt/gc2719/halllab/bin/bwa
	# Description	: alignment via Burrows-Wheeler transformation
	# Version		: 0.7.8-r455
	# Author			: Heng Li <lh3@sanger.ac.uk>
	# Usage			: speedseq realign $BAM
	# Command		: index	index sequences in the FASTA format
	#          		  mem	BWA-MEM algorithm
# 2. SAMBLASTER
	# Program 		: samblaster
	# Local Path			: /gscmnt/gc2719/halllab/bin/samblaster
	# Version		: 0.1.21
	# Author			: Greg Faust (gf4ea@virginia.edu)
	# Summary		: Tool to mark duplicates and optionally output split reads and/or discordant pairs.
	# Usage			: e.g. 1. bwa mem index samp.r1.fq samp.r2.fq | samblaster [-e] [-d samp.disc.sam] [-s samp.split.sam] | samtools view -Sb - > samp.out.bam
	# 	     			  		 2. samtools view -h samp.bam | samblaster [-a] [-e] [-d samp.disc.sam] [-s samp.split.sam] [-u samp.umc.fasta] -o /dev/null
# 3. SAMBAMBA 
	# Program 		: sambamba
	# Local Path			: /gscmnt/gc2719/halllab/bin/sambamba
	# Version		: v0.5.4
	# Author			: Artem Sarasov
	# Summary		: Faster and parallel implementation to work with SAM and BAM files
	# Usage			: available commands : 'view', 'index', 'merge', 'sort', 'flagstat', 'slice', 'markdup', 'depth', 'mpileup' 
# 4  NUMPY and Scipy
	# Package 		: Numpy version 1.8.1; Scipy version 0.14
	# Python fundamental packages for scientific computing with Python
# 5  pysam 
		# Package	: pysam version 0.8.0+ 
		# Python package for working with alignment files in SAM/BAM format 	
# 6. ROOT 
	# Program 		: ROOT
	# Author			: CERN
	# Summary		: ROOT is a framework for data processing: saving, mining, accessing data, and enables publishing results in graphs
# 7. Variant Effect Predictor 
	# Program 		: VEP
	# Author			: ENSEMBL
	# Summary		: VEP determines the effect of variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions
# 8. VAWK:
	# Program 		: vawk
	# Local Path			: /gscmnt/gc2719/halllab/bin/vawk
	# Version		: 0.0.2
	# Author			: Colby Chiang
	# Summary		: An awk-like VCF parser
	# Usage			: usage: vawk [-h] [-v VAR] [-c INFO_COL] [--header] [--debug] cmd [vcf]
	
# 9. MBUFFER:
	# Program		: mbuffer
	# Summary		: mbuffer is a tool for buffering data streams over TCP based network targets, paralleling stream, etc
	# Usage 			: mbuffer -b<num> -s<size> -m<size> -i<file> -o<file>
# 10. GNU Parallel
	# Program		: parallel
	# Summary 		: a shell tool for executing jobs in parallel
 

#NOTE: All dependencies need to be installed before speedseq is installed. 
#Detailed instructions for speedseq installation and dependencies at: https://github.com/hall-lab/speedseq
# (2) construct the alignment command
for SAMPLE in `cat $SAMPLEMAP | awk '{ print $1 }'`
do
    ORIG_BAM=`cat $SAMPLEMAP | awk -v SAMPLE=$SAMPLE '{FS="\t"; if ($1==SAMPLE) print $2 }'`
    echo "mkdir -p $WORKDIR/notes/$BATCH/log && \
        mkdir -p $WORKDIR/$SAMPLE && \
        bomb -q hall-lab -m 48 -t 8 -J $SAMPLE \
        -o $WORKDIR/notes/$BATCH/log/$SAMPLE.align.%J.log \
        -e $WORKDIR/notes/$BATCH/log/$SAMPLE.align.%J.log \
         \"speedseq realign \
            -o $WORKDIR/$SAMPLE/$SAMPLE \
            -T $WORKDIR/$SAMPLE/temp \
            -M 8 \
            -t 8 \
            -v \
            /gscmnt/gc2719/halllab/genomes/human/GRCh37/hs37_ebv/hs37_ebv.fasta \
            $ORIG_BAM\""
done > $WORKDIR/notes/realign_command.sh


# ----------------------------------------
# 3. run alignment
# ----------------------------------------
cat $WORKDIR/notes/realign_command.sh | bash

# ----------------------------------------
# 4. Flagstat the realigned BAM files
# ----------------------------------------
for SAMPLE in `cat $BATCH | cut -f 1`
do
    bomb \
        -m 1 \
        -J $SAMPLE.flag \
        -o $WORKDIR/$SAMPLE/log/$SAMPLE.realign.flagstat.%J.log \
        -e $WORKDIR/$SAMPLE/log/$SAMPLE.realign.flagstat.%J.log \
        "sambamba flagstat $WORKDIR/$SAMPLE/$SAMPLE.bam > $WORKDIR/$SAMPLE/$SAMPLE.bam.flagstat"
done

#SET DIRECTORY HERE:
DIR="/gscmnt/gc2802/halllab/users/abadve/projects/illumina_platinum/"
mkdir $DIR/notes/

# create a delimited batch file with sample and bam paths 
# e.g. NA12877	/gscmnt/gc2802/halllab/ceph1463_realign_021815/NA12877/NA12877.bam

for SAMPLE in `cat $DIR/notes/batch.txt | cut -f 1`
do
    mkdir -p $DIR/lumpy/$SAMPLE/log
done

# ---------------------------------------------------------------------
# 5a. Run LUMPY through Speedseq SV to call structural variants
# ---------------------------------------------------------------------
# ----------------------------------------
# speedseq
# ----------------------------------------
# Program: speedseq
# Path : /gscmnt/gc2719/halllab/bin/speedseq
# Version: 0.0.3
# Author: Colby Chiang
# usage:   speedseq <command> [options]
# command: align    align FASTQ files with BWA-MEM 
#          var      call SNV and indel variants with FreeBayes
#          somatic  call somatic SNV and indel variants in a tumor/normal pair with FreeBayes
#          sv       call SVs with LUMPY
#          realign  realign from a coordinate sorted BAM file 
#
# LUMPY    
	# Program   	: lumpy
	# Path			: /gscmnt/gc2719/halllab/bin/lumpy
	# Version		: v 0.2.11
	# Summary		: Find structural variations in various signals.
	# Author			: Ryan Layer
	# 	Usage	   	: speedseq sv [options]

while read SAMPLE BAM PROJECT 
do 
SPL="${BAM%.*}.splitters.bam"
DCD="${BAM%.*}.discordants.bam"
bomb \
-g /abadve/lumpy \
-m 30 \
-J $SAMPLE.lumpy \
-o $DIR/lumpy/$SAMPLE/log/$SAMPLE.lumpy.log \
-e $DIR/lumpy/$SAMPLE/log/$SAMPLE.lumpy.log \
"speedseq sv \
-o $DIR/lumpy/$SAMPLE/$SAMPLE \
-T $DIR/lumpy/$SAMPLE/temp \
-R /gscmnt/gc2719/halllab/genomes/human/GRCh37/hs37_ebv/hs37_ebv.fasta \
-B $BAM \
-S $SPL \
-D $DCD \
-v \
-x /gscmnt/gc2719/halllab/src/speedseq/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed \
-d \
-P \
-g \
-k"
done < $DIR/notes/batch.txt

# ------------------------------------------------------------
# 5b. Count SVs after lumpy
# ------------------------------------------------------------
#After running lumpy to detect structural variants, 
#run following command to count various SV types
#tandem duplications, deletions,inversions and break-ends

for SAMPLE in `cat $DIR/notes/batch.txt | cut -f 1` 
do
zcat $DIR/lumpy/$SAMPLE/$SAMPLE.sv.vcf.gz \
| vawk -v S=$SAMPLE 'BEGIN {DEL=0; DUP=0; INV=0; BND=0} { if (I$SVTYPE=="DEL") DEL+=1; \
else if (I$SVTYPE=="DUP") DUP+=1; else if (I$SVTYPE=="INV") INV+=1; \
else if (I$SVTYPE=="BND" && ! I$SECONDARY) BND+=1 } END { print S,DEL,DUP,INV,BND }'
done >> $DIR/notes/svtype_counts.txt

# ---------------------------------------------
# 6a. Genotyping the SV Callset using genotype
# ---------------------------------------------
# Lumpy outputs compressed vcf output files which need to be extracted 
# for running genotype.
# Following command extracts bunch of compressed(.gz) vcf files into raw vcf format
for SAMPLE in `cat $DIR/notes/batch.txt | cut -f 1`
do
    echo $SAMPLE
    zcat $DIR/lumpy/$SAMPLE/$SAMPLE.sv.vcf.gz > $DIR/lumpy/$SAMPLE/$SAMPLE.sv.vcf
done


# -----------------------------------------------------------
# 6b. Using lsort we want to sort and concatenate VCF files
# ----------------------------------------------------------

# ----------------------------------------
# lsort  
# ----------------------------------------
# Program: lsort
# Author: Ryan Layer and Ira Hall
# Path: /gscmnt/gc2719/halllab/bin/lsort
# Version: 0.01 
# Description: sort N VCF files into a single file
# Usage: %prog <VCF file 1> <VCF file 2> ... <VCF file N>
# 
echo -n svtools lsort > $DIR/sort_cmd.sh
for SAMPLE in `cat $DIR/notes/batch.txt | cut -f 1`
do
    echo -ne " \\\\\n\t$DIR/lumpy/$SAMPLE/$SAMPLE.sv.vcf"
done >> $DIR/sort_cmd.sh
bomb -m 25 -J lsort "bash $DIR/sort_cmd.sh | bgzip -c > $DIR/sorted.sv.vcf.gz"

# -----------------------------------------------------------
# 6c. Collapse the variants into merged VCF
# ----------------------------------------------------------

# ----------------------------------------
# lmerge  
# ----------------------------------------
# Program: lmerge
# Author: Ryan Layer and Ira Hall
# Path: /gscmnt/gc2719/halllab/bin/lmerge
# Version: ira_7 
# Description:merge lumpy calls.
# Usage: lmerge -i <file>
# Dependencies on other software:
# 1. l_bp.py: Implementation to decide which variants can be merged 

bomb -m 20 -J lmerge.$SLOP \
    "zcat $DIR/sorted.sv.vcf.gz \
        |svtools lmerge -i /dev/stdin -f 20 \
        | bedtools sort -header \
        | bgzip -c \
        > $DIR/merged.sv.vcf.gz"

# ----------------------------------------
# 6d. Genotype merged VCF with genotype
# ----------------------------------------
# To speed things up, generate a separate VCF for each sample and
# join them afterwards.

# ----------------------------------------
# Genotype (earlier called SVTYPER)
# ----------------------------------------
# Program: genotype
# Author: Colby Chiang 
# Path: /gscmnt/gc2719/halllab/bin/genotype
# Version: 0.0.2 
# Description: Compute genotype of structural variants based on breakpoint depth
# Usage:  genotype [-h] -B BAM -S SPLIT_BAM [-i INPUT_VCF] [-o OUTPUT_VCF] [-f SPLFLANK] [-F DISCFLANK] [--split_weight SPLIT_WEIGHT][--disc_weight DISC_WEIGHT] [-n NUM_SAMP] [--debug]

#Dependencies:
#  pysam 
		# Package	: pysam version 0.8.0+ 
		# Python package for working with alignment files in SAM/BAM format 	

mkdir -p gt
while read SAMPLE BAM 
do 
    echo $SAMPLE
    SPL=${BAM%.*}.splitters.bam
    bomb -m 18 -J $SAMPLE.gt -o log/$SAMPLE.gt.%J.log -e log/$SAMPLE.gt.%J.log \
         "zcat merged.sv.vcf.gz \
            | vawk --header '{  \$6=\".\"; print }' \
            | svtools genotype \
                -B $BAM \
                -S $SPL \
            | sed 's/PR...=[0-9\.e,-]*\(;\)\{0,1\}\(\t\)\{0,1\}/\2/g' - \
            > gt/$SAMPLE.vcf"
done < $DIR/notes/batch.txt

# ----------------------------------------
# 6e. paste the samples into a single VCF
# ----------------------------------------
# ----------------------------------------
# Vcfpaste  
# ----------------------------------------
# Program: vcfpaste
# Author: Colby Chiang / Abhijit Badve
# Path: /gscmnt/gc2719/halllab/bin/vcfpaste
# Version: 0.0.2 
# Description:Paste VCFs from multiple samples
# Usage: vcfpaste [-h] [-m MASTER] [vcf [vcf ...]]

bomb -m 20 -J paste.gt -o log/paste.gt.%J.log -e log/paste.gt.%J.log \
    "zcat merged.sv.vcf.gz \
        | vcf_group_multiline.py \
        | svtools vcfpaste \
            -m - \
            gt/*.vcf \
        | bedtools sort -header \
        > merged.sv.gt.vcf"



# ---------------------------------------------------------
# 7. Annotate read-depth of VCF variants
# add copy number information to the VCF using copynumber
# ---------------------------------------------------------

# ----------------------------------------
# CopyNumber  
# ----------------------------------------
# Program: copynumber
# Author: Colby Chiang
# Path: /gscmnt/gc2719/halllab/bin/copynumber
# Version: 0.0.1 
# Description: Compute genotype of structural variants based on breakpoint depth
# Usage: copynumber [-h] [-v INPUT_VCF] [-r ROOT] [-w WINDOW] [-s SAMPLE] [--cnvnator CNVNATOR] [-o OUTPUT_VCF] [--debug]
# Calls CNVNATOR-MULTI inside to get the copy number information 
 # To speed things up, generate a separate VCF for each sample and
# join them afterwards.
mkdir -p cn
#ROOT libraries path
source /gsc/pkg/root/root/bin/thisroot.sh
VCF=merged.sv.gt.vcf
create_coordinates.py -i $VCF -o coordinates
while read SAMPLE BAM PROJECT 
do
    ROOT=$DIR/lumpy/$SAMPLE/temp/cnvnator-temp/$SAMPLE.bam.hist.root
    SM=`sambamba view -H $BAM | grep -m 1 "^@RG" | awk '{ for (i=1;i<=NF;++i) { if ($i~"^SM:") SM=$i; gsub("^SM:","",SM); } print SM }'`
    bomb -q hall-lab -m 6 -J $SAMPLE.cn -o log/$SAMPLE.cn.%J.log -e log/$SAMPLE.cn.%J.log \
         "svtools copynumber \
             --cnvnator cnvnator-multi \
             -s $SM \
             -w 100 \
             -r $ROOT \
             -c coordinates \
             -v $VCF \
          > cn/$SAMPLE.vcf"
done < $DIR/notes/batch.txt

# paste the samples into a single VCF
bomb -m 4 -J paste.cn -o log/paste.cn.%J.log -e log/paste.cn.%J.log \
    "svtools vcfpaste \
        -m merged.sv.gt.vcf \
        cn/*.vcf \
        | bgzip -c \
        > merged.sv.gt.cn.vcf.gz"

# -----------------------------------------------------
# 8. To accurately count SVs we need to further prune 
# sv callsets to find overlapping variants within a vcf
# Criteria to pruning variants is minimum overlapping 
# distance(Default 50 ) and comparing eval parameter 
# e.g.(allele frequencies)
# -----------------------------------------------------
# To run prune we need to convert vcf file to bedpe 
# Also if evaluating parameter is allele frequency we need to add AF using afreq 
# -s option in prune specifies if the file is sorted to minimize stack usage

# ----------------------------------------
# Prune 
# ----------------------------------------
# Program: prune
# Author: Abhijit Badve 
# Path: /gscmnt/gc2719/halllab/bin/prune
# Version: 0.0.1 
# Description:  cluster a BEDPE file by position based on their allele frequency
# Usage: prune [-h] [-d INT] [-e string] [-s] [-o OUTPUT] [input]

# # ----------------------------------------
# Afreq  
# ----------------------------------------
# Program: afreq
# Author: Colby Chiang
# Path: /gscmnt/gc2719/halllab/bin/afreq
# Version: 0.0.1 
# Description: Add allele frequency information to a VCF file
# Usage: afreq [-h] [vcf]

#command
zcat merged.sv.gt.cn.vcf.gz | svtools vcftobedpe \
	| svtools afreq \
	| svtools prune -d 100 -e "AF" -s |  svtools bedpetovcf -o merged.pruned.sv.gt.cn.vcf.gz

# ---------------------------------------------------
# 9. Compare SV counts for pre- and post-merged VCFs
# ---------------------------------------------------

# count up the pre-merged
mkdir -p pre-merged_sv_count
while read SAMPLE BAM 
do 
    SM=`sambamba view -H $BAM | grep -m 1 "^@RG" | awk '{ for (i=1;i<=NF;++i) { if ($i~"^SM:") SM=$i; gsub("^SM:","",SM); } print SM }'`
	 bomb  -m 2 -J $SAMPLE.svcount  -o log/%J.log -e log/%J.log \
        "svcounts.sh \
            pre-merged_gt/$SAMPLE.sv.gt.vcf.gz $SM 0 \
            > pre-merged_sv_count/$SAMPLE.count.q0.txt"
done < $DIR/notes/batch.txt

# consolidate into a single file
cat pre-merged_sv_count/*.count.q0.txt \
    | awk '{ print $1,$2,$4+$5 }' OFS="\t" \
    | groupBy -g 1 -c 3 -o collapse \
    | tr ',' '\t' \
    > pre-merged.sv_counts.q0.txt

# count up the post-merged
mkdir -p sv_count
while read SAMPLE BAM PROJECT 
do 
	SM=`sambamba view -H $BAM | grep -m 1 "^@RG" | awk '{ for (i=1;i<=NF;++i) { if ($i~"^SM:") SM=$i; gsub("^SM:","",SM); } print SM }'`
	bomb -m 2 -J $SAMPLE.svcount -o log/%J.log -e log/%J.log \
	    "sv_counts.sh \
	        gt/$SAMPLE.vcf $SM 0 \
	        > sv_count/$SAMPLE.count.q0.txt"
done < $DIR/notes/batch.txt

# consolidate into a single file
cat sv_count/*.count.q100.txt \
    | awk '{ print $1,$2,$4+$5 }' OFS="\t" \
    | groupBy -g 1 -c 3 -o collapse \
    | tr ',' '\t' \
    > post-merged.sv_counts.q100.txt
	
# ----------------------------------------------------------
# 10. VarLookup: To compare multiple subsets and discover
#variants common(overlapping at min distance) between them. 
# ----------------------------------------------------------
# ----------------------------------------
# varlookup 
# ----------------------------------------
# Program: varlookup
# Author: Abhijit Badve 
# Path: /gscmnt/gc2719/halllab/bin/varlookup
# Version: 0.0.1 
# Description:Look for variants common between two bedpe files
# Usage: varlookup [-h] [-d INT] [-a FILE] [-b FILE] [-c string] [-o OUTPUT]

#command
zcat merged.pruned.cohort1.sv.gt.cn.bedpe | svtools varlookup -a stdin -b merged.pruned.cohort2.sv.gt.cn.bedpe \
> merged.pruned.cohort1_2.sv.gt.cn.bedpe

# ----------------------------------------------------------
# 11. Classify: Classify structural variants based on an 
#annotation 
# ----------------------------------------------------------
# ----------------------------------------
# Classify 
# ----------------------------------------
# author: Abhijit Badve
# version: $Revision: 0.0.1 $
# description: classify structural variants
#
# optional arguments:
#   -h, --help            show this help message and exit
#   -t String, --tSet String
#                         high quality deletions & duplications training dataset[vcf]/[stdin]
#   -i VCF, --input VCF   test vcf input for applying a model
#   -o VCF, --output VCF  vcf output [stdout]
#   --debug               debugging verbosity
#
#Step 1:
#Create a training VCF
#fixed set of known, high-quality variants that then typed for each set of individuals to train the model
#Step 2:
#command
cat training.data.vcf | classify -t - -i test.vcf  -o out.vcf
#Note: training vcf should have all the individual samples as in test vcf

# ----------------------------------------------------
# 9. Plot the SV counts in R before and after merging
# ----------------------------------------------------
#Following part of the tutorial is run in R environment
mkdir -p plots
# R code:
# 
col.list <- c('indianred3', 'dodgerblue3', 'goldenrod1', 'gray')

#Read svcounts before merging
before <- read.table('pre-merged.sv_counts.q0.txt')
colnames(before)<-c("SAMPLE","DEL","DUP","INV","BND","PROJECT")
customlevels.before<-as.matrix(unique(before$PROJECT))
before.ordered <- before[order(factor(before$PROJECT,levels = customlevels.before)),]

#Read svcounts after merging
after <- read.table('post-merged.sv_counts.q0.txt')
colnames(after)<-c("SAMPLE","DEL","DUP","INV","BND","PROJECT")
customlevels.after<-as.matrix(unique(after$PROJECT))
after.ordered <- after[order(factor(after$PROJECT,levels = customlevels.after)),]

chngpnt<-data.frame()
unname(chngpnt)

unname(before.ordered)
for (i in 2:(nrow(before.ordered))) {  
   if (before.ordered[i,"PROJECT"]!=before.ordered[(i-1),"PROJECT"]) {  
		chngpnt<-rbind(chngpnt,i-1)
     }  
 } 
chngpnt<-rbind(chngpnt,i)
row.names(chngpnt)<-customlevels.before

chngpnt<-data.frame()
unname(chngpnt)
unname(after.ordered)
for (i in 2:(nrow(after.ordered))) {  
   if (after.ordered[i,"PROJECT"]!=after.ordered[(i-1),"PROJECT"]) {  
		chngpnt<-rbind(chngpnt,i-1)
     }  
 } 
chngpnt<-rbind(chngpnt,i)
row.names(chngpnt)<-customlevels.after



# before

pdf('plots/pre-merged_sv_counts.q0.pdf', height=10, width=20)
par(lty=0)
barplot(t(data.matrix(before.ordered[,2:(ncol(before.ordered)-1)])), col=col.list, las=3, cex.names=0.6, main="before merging\nnumber of SVs per sample by type with no score-filter", xaxt='n',ylab='Number of SVs', ylim=c(0,16000), xaxt='n', space=0)
#barplot(t(data.matrix(before[,2:ncol(before)])), col=col.list, las=3, cex.names=0.6, main="Before merging\nnumber of SVs per sample by type", ylab='Number of SVs', ylim=c(0,9000), xlab='Samples')
legend('topleft', c('Deletions', 'Tandem duplications', 'Inversions', 'Unknown'), fill=col.list, bty='n')
a=0
par(lty=4)
for(i in 1:length(customlevels.before)) {
    b=chngpnt[customlevels.before[i],]
    abline(v=b,col="grey")
    print(paste((a+b)/2,b,sep=" "))
  	 text(x = (a+b)/2,y = 6500, labels=customlevels.before[i], srt=90,cex = 0.9,col="blue",lty=4)
    a=b
   
}
dev.off()
# after
pdf('plots/post-merged_sv_count.q0.pdf', height=10, width=20)
par(lty=0)
barplot(t(data.matrix(after.ordered[,2:(ncol(after.ordered)-1)])), col=col.list, las=3, cex.names=0.6, main="after merging\nnumber of SVs per sample by type with no score-filter",  xaxt='n', ylab='Number of SVs', ylim=c(0,16000), space=0)
legend('topleft', c('Deletions', 'Tandem duplications', 'Inversions', 'Unknown'), fill=col.list, bty='n')
a=0
par(lty=4)
for(i in 1:length(customlevels.after)) {
    b=chngpnt[customlevels.after[i],]
    abline(v=b,col="grey")
    print(paste((a+b)/2,b,sep=" "))
    text(x = (a+b)/2,y = 6500, labels=customlevels.after[i], srt=90,cex = 0.9,col="blue",lty=4)
    a=b
   
}
dev.off()
# end R code


#PROGRAMS USED IN THIS TUTORIAL: 
# ----------------------------------------
# 1. speedseq
# ----------------------------------------
# Program: speedseq
# Path : /gscmnt/gc2719/halllab/bin/speedseq
# Version: 0.0.3
# Author: Colby Chiang
# usage:   speedseq <command> [options]
# command: align    align FASTQ files with BWA-MEM 
#          var      call SNV and indel variants with FreeBayes
#          somatic  call somatic SNV and indel variants in a tumor/normal pair with FreeBayes
#          sv       call SVs with LUMPY
#          realign  realign from a coordinate sorted BAM file 
#
# Dependencies on other software:
# 1. BWA: speedseq uses bwa mem to   
	# Program   : bwa mem:
	# Path			: /gscmnt/gc2719/halllab/bin/bwa
	# Description	: alignment via Burrows-Wheeler transformation
	# Version		: 0.7.8-r455
	# Author			: Heng Li <lh3@sanger.ac.uk>
	# Usage			: speedseq realign $BAM
	# Command		: index	index sequences in the FASTA format
	#          		  mem	BWA-MEM algorithm
# 2. LUMPY    
	# Program   	: lumpy
	# Path			: /gscmnt/gc2719/halllab/bin/lumpy
	# Version		: v 0.2.11
	# Summary		: Find structural variations in various signals.
	# Author			: Ryan Layer
	# 	Usage	   	: speedseq sv [options]
# 3. freebayes
	# Program		: freebayes
	# Path			: /gscmnt/gc2719/halllab/bin/freebayes
	# version		: v0.9.21-7-g7dd41db
	# Summary		: Bayesian haplotype-based polymorphism discovery.
	# author			: Erik Garrison, Marth Lab, Boston College, 2010-2014
	# usage			: [speedseq var [REFERENCE] [OPTIONS] [BAM FILES] >[OUTPUT]
# 4. SAMBLASTER
	# Program 		: samblaster
	# Path			: /gscmnt/gc2719/halllab/bin/samblaster
	# Version		: 0.1.21
	# Author			: Greg Faust (gf4ea@virginia.edu)
	# Summary		: Tool to mark duplicates and optionally output split reads and/or discordant pairs.
	# Usage			: e.g. 1. bwa mem index samp.r1.fq samp.r2.fq | samblaster [-e] [-d samp.disc.sam] [-s samp.split.sam] | samtools view -Sb - > samp.out.bam
	# 	     			  		 2. samtools view -h samp.bam | samblaster [-a] [-e] [-d samp.disc.sam] [-s samp.split.sam] [-u samp.umc.fasta] -o /dev/null
# 5. SAMBAMBA 
	# Program 		: sambamba
	# Path			: /gscmnt/gc2719/halllab/bin/sambamba
	# Version		: v0.5.4
	# Author			: Artem Sarasov
	# Summary		: Faster and parallel implementation to work with SAM and BAM files
	# Usage			: available commands : 'view', 'index', 'merge', 'sort', 'flagstat', 'slice', 'markdup', 'depth', 'mpileup'

# ----------------------------------------
# 2. genotype (SVTYPER)
# ----------------------------------------
# Program: genotype
# Author: Colby Chiang 
# Path: /gscmnt/gc2719/halllab/bin/genotype
# Version: 0.0.2 
# Description: Compute genotype of structural variants based on breakpoint depth
# Usage:  genotype [-h] -B BAM -S SPLIT_BAM [-i INPUT_VCF] [-o OUTPUT_VCF] [-f SPLFLANK] [-F DISCFLANK] [--split_weight SPLIT_WEIGHT][--disc_weight DISC_WEIGHT] [-n NUM_SAMP] [--debug]

# ----------------------------------------
# 3. Prune 
# ----------------------------------------
# Program: prune
# Author: Abhijit Badve 
# Path: /gscmnt/gc2719/halllab/bin/prune
# Version: 0.0.1 
# Description:  cluster a BEDPE file by position based on their allele frequency
# Usage: prune [-h] [-d INT] [-e string] [-s] [-o OUTPUT] [input]

# ----------------------------------------
# 4. Varlookup 
# ----------------------------------------
# Program: varlookup
# Author: Abhijit Badve 
# Path: /gscmnt/gc2719/halllab/bin/varlookup
# Version: 0.0.1 
# Description:Look for variants common between two bedpe files
# Usage: varlookup [-h] [-d INT] [-a FILE] [-b FILE] [-c string] [-o OUTPUT]

# ----------------------------------------
# 5. Classify 
# ----------------------------------------
# Program: classify
# Author: Colby Chiang
# Path: /gscmnt/gc2719/halllab/bin/classify
# Version: 0.0.2 
# Description:classify structural variants
# Usage: classify [-h] [-i VCF] [-g FILE] [-e FILE] [-a BED] [-f FLOAT] [-s FLOAT] [-r FLOAT]

# ----------------------------------------
# 6. Vcfpaste  
# ----------------------------------------
# Program: vcfpaste
# Author: Colby Chiang / Abhijit Badve
# Path: /gscmnt/gc2719/halllab/bin/vcfpaste
# Version: 0.0.2 
# Description:Paste VCFs from multiple samples
# Usage: vcfpaste [-h] [-m MASTER] [vcf [vcf ...]]

# ----------------------------------------
# 7. CopyNumber  
# ----------------------------------------
# Program: copynumber
# Author: Colby Chiang
# Path: /gscmnt/gc2719/halllab/bin/copynumber
# Version: 0.0.1 
# Description: Compute genotype of structural variants based on breakpoint depth
# Usage: copynumber [-h] [-v INPUT_VCF] [-r ROOT] [-w WINDOW] [-s SAMPLE] [--cnvnator CNVNATOR] [-o OUTPUT_VCF] [--debug]

# ----------------------------------------
# 8. Afreq  
# ----------------------------------------
# Program: afreq
# Author: Colby Chiang
# Path: /gscmnt/gc2719/halllab/bin/afreq
# Version: 0.0.1 
# Description: Add allele frequency information to a VCF file
# Usage: afreq [-h] [vcf]

# ----------------------------------------
# 9. Lmerge  
# ----------------------------------------
# Program: lmerge
# Author: Ryan Layer and Ira Hall
# Path: /gscmnt/gc2719/halllab/bin/lmerge
# Version: ira_7 
# Description:merge lumpy calls.
# Usage: lmerge -i <file>
# Dependencies on other software:
# 1. l_bp.py: Implementation to decide which variants can be merged 

# ----------------------------------------
# 10. Lsort  
# ----------------------------------------
# Program: lsort
# Author: Ryan Layer and Ira Hall
# Path: /gscmnt/gc2719/halllab/bin/lsort
# Version: 0.01 
# Description: sort N VCF files into a single file
# Usage: %prog <VCF file 1> <VCF file 2> ... <VCF file N>

# ----------------------------------------
# 10. Bomb  
# ----------------------------------------
# Program: bomb
# Author: Colby Chiang
# Path: /gscmnt/gc2719/halllab/bin/bomb
# Version: 0.0.1 
# Description: Constructor for bsub jobs
# Usage: [-h] [-J NAME] [-g GROUP] [-q QUEUE] [-t INT] [-m INT] [-n INT][-x HOST] [-i HOST] [-o OUT] [-e ERR] [-E EMAIL] [-d] "COMMAND"

