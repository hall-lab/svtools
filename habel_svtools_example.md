#Example Analysis using svtools
This tutorial will help you begin to explore the use of svtools to analyze an SV vcf.  It will help you to satisfy the computing environment requirements, gather the required genomic data, and try basic analysis using svtools.
This tutorial includes example commands that you can alter to refer to your sample names.
```
:note A shell script that implements this tutorial for a small set of samples has also been included
```
[lumpy_pipeline.sh](https://github.com/jeldred/svtools/edit/install_documentation/lumpy_pipeline.sh)
```
Creating a sample.map file and running lumpy_pipeline.sh is another way to investigate the usage 
of svtools for creating a callset.
```
##Satisfy computing environment requirements
1. Install svtools (and Python)
2. Acquire vawk - copy vawk into demo directory or  (/usr/bin/vawk)?
3. Install CNVnator

###Install svtools (and Python)
[Installation instructions](https://github.com/jeldred/svtools/blob/install_documentation/INSTALL.md) have been provided in the INSTALL.md and DEVELOPER.MD of this repo.
###Acquire vawk
clone from the repo? (how should we recommend people get this in their path and test to make sure it is falling through to the right gawk) [vawk](https://github.com/cc2qe/vawk) is an awk-like VCF parser that is used to manipulate vcf files in this tutorial.
###Install CNVnator(?)
for use with cnvnator, can we change this to link them out to generally prepare to run cnvnator?  and then maybe just note the need to source thisroot.sh below somewhere [CNVnator](https://github.com/abyzovlab/CNVnator/releases)
```
note: thisroot.sh is required for the Copy Number Annotation step in this tutorial.
```

##Gather genomic data and generate needed helper files

2. Get or Create SpeedSeq/lumpy SV vcf files
3. Get or Create SpeedSeq aligned bams and spliiter files
4. Refernce FASTA 
5. cn.list file

### Get or Create SpeedSeq/lumpy SV vcf files
### Get or Create aligned bams and splitter files
### Reference FASTA
### cn.list file
The cn.list file has a single column that contains the path to the vcf files output in the Copy Number Annotation step of this tutorial.

##Use svtools to "process" data FIXME
1. Use vawk to remove 'REF' variants from SpeedSeq SV vcf
2. Use svtools lsort to combine and sort variants from multiple samples
3. Use svtools lmerge to merge variants deemed to be identical in the sorted vcf
4. Remove EBV (Epstein-Barr Virus) variants
5. Use svtools genotype to force genotypes for variant positions discovered in other samples.
6. Copy Number annotation - ROOT libraries path 
    1. prepare environemnt for cnvnator
    2. make uncompressed copy 
    3. make coordinate file
    4. Annotate variants with copy number from cnvnator using FIXME
7. Use svtools vcfpaste to FIXME
8. Prune pipeline 
9. Training (exclude?)
10. Classification (exclude?)

###Use vawk to remove 'REF' variants from SpeedSeq SV vcf
This step will remove variants that have been detected but then determined to be homozygous reference.
```
  zcat SAMPLE1.sv.vcf.gz \
  | vawk --header '{if(S$*$GT!="0/0" && S$*$GT!="./.") print $0}' \
  > SAMPLE1.sv.non_ref.vcf
```

### Use svtools lsort to combine and sort variants from multiple samples
```
svtools lsort SAMPLE1.sv.non_ref.vcf SAMPLE2.sv.non_ref.vcf SAMPLE3.sv.non_ref.vcf \
| bgzip -c > sorted.vcf.gz
```
```
:note svtools lsort will remove variants with the SECONDARY tag in the INFO field.
This will cause the sorted vcf to have fewer variant lines than the input.
```
### Use svtools lmerge to merge variants deemed to be identical in the sorted vcf
```
zcat sorted.vcf.gz \
  | svtools lmerge -i /dev/stdin --product -f 20 \
  | bgzip -c > merged.vcf.gz "

```
```
:note svtools lmerge will return variant lines for SECONDARY break ends in addition to merging variants.
This will sometimes cause the merged vcf to have more variant lines than the input.
```
### Remove variants detected by alignment to the EBV (Epstein-Barr Virus) contig
```
zcat merged.vcf.gz \
| $VAWK --header '{if($0 !~ /NC_007605/) print $0}' \
| bgzip -c > merged.no_EBV.vcf.gz
```

### Use svtools genotype to force genotypes for variant positions discovered in other samples
svtools genotype will calculate a genotype for each sample at the variant positions in the merged.no_EBV.vcf.gz file.
It requires the aligned bam and splitters file for each sample.
```
    "zcat merged.no_EBV.vcf.gz \
     | vawk --header '{  \$6=\".\"; print }' \
     | svtools genotype \
       -B SAMPLE1.bam \
       -S SAMPLE1.splitters.bam \
     | sed 's/PR...=[0-9\.e,-]*\(;\)\{0,1\}\(\t\)\{0,1\}/\2/g' - \
     > SAMPLE1.vcf"
```

###Copy Number annotation - ROOT libraries path 
####prepare directory structure
mkdir -p cn/logs/MISC
mkdir -p cn/MISC
#### prepare environemnt for cnvnator
```
source /gsc/pkg/root/root/bin/thisroot.sh
```
#### make uncompressed copy 
```
zcat merged.no_EBV.vcf.gz > merged.no_EBV.vcf
VCF=merged.no_EBV.vcf
```
#### make coordinate file
```
create_coordinates -i $VCF -o coordinates

while read SAMPLE BAM
do
  base=`basename $BAM .bam`
  ROOT=/gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/$SAMPLE/temp/cnvnator-temp/$base.bam.hist.root
  # cnvnator files were generated by speedseq that I did not run as part of this demo

  #SM=`sambamba view -H $BAM | grep -m 1 "^@RG" | awk '{ for (i=1;i<=NF;++i) { if ($i~"^SM:") SM=$i; gsub("^SM:","",SM); } print SM }'`

echo SAMPLE $SAMPLE
#echo SM $SM

  bsub -M 4000000  -R 'select[mem>4000] rusage[mem=4000]' -u jeldred@genome.wustl.edubsub -q long -J $SAMPLE.cn -o cn/logs/$SAMPLE.cn.%J.log -e cn/logs/$SAMPLE.cn.%J.log \
     "svtools copynumber \
         --cnvnator /gscmnt/gc2719/halllab/bin/cnvnator-multi \
         -s $SAMPLE \
         -w 100 \
         -r $ROOT \
         -c coordinates \
         -v gt/$SAMPLE.vcf \
      > cn/$SAMPLE.vcf"
done < sample.map
```

### VCF paste
```
VCF=merged.no_EBV.vcf
bsub -M 4000000  -R 'select[mem>4000] rusage[mem=4000]' -u jeldred@genome.wustl.edubsub -q long -J paste.cn -o paste.cn.%J.log -e paste.cn.%J.log \
    "svtools vcfpaste \
        -m $VCF \
        -f cn.list \
        -q \
        | bgzip -c \
        > merged.sv.gt.cn.vcf.gz"
```

### Prune pipeline 
```
bsub -q long -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' "zcat merged.sv.gt.cn.vcf.gz \
| svtools afreq \
| svtools vcftobedpe \
| svtools bedpesort \
| svtools prune -s -d 100 -e \"AF\" \
| svtools bedpetovcf \
| bgzip -c > merged.sv.new_pruned.vcf.gz"
```

### Training
```
zcat merged.sv.new_pruned.vcf.gz \
 | svtools vcftobedpe  \
 | svtools varlookup -a stdin -b /gscmnt/gc2802/halllab/sv_aggregate/reclass/finmetseq.training_vars.bedpe.gz -c FINMETSEQ_HQ -d 50 \
 | svtools bedpetovcf \
 | $VAWK --header '{if(I$FINMETSEQ_HQ_AF>0) print $0}' \
 | bgzip -c > training.vars.vcf.gz
```

## Reclassify THIS SECTION WELL OUT OF DATE FIXME
```
zcat merged.sv.new_pruned.vcf.gz \
 | python /gscmnt/gc2802/halllab/sv_aggregate/dev/svtools/svtools/reclass_combined.py -g /gscmnt/gc2802/halllab/sv_aggregate/ceph_ped/ceph.sex.txt  -t <(zcat training.vars.vcf.gz)  -a /gscmnt/gc2719/halllab/users/cchiang/projects/g#tex/annotations/repeatMasker.recent.lt200millidiv.b37.sorted.bed.gz  -d class.diags.0313.txt  | bgzip -c > reclass.0313.all.vcf.gz



#zcat reclass.0313.all.vcf.gz \
#| vawk --header '{
#  split(I$STRANDS,x,",");
#  split(x[1],y,":");
#  split(x[2],z,":");
#  if (I$SVTYPE=="DEL" || I$SVTYPE=="DUP" || I$SVTYPE=="MEI"){
#   $7="PASS"; print $0;
#  }  else if ( I$SVTYPE=="INV" && $8>=100 && (I$SR/I$SU)>=0.1 && (I$PE/I$SU)>=0.1 && (y[2]/I$SU)>0.1 && (z[2]/I$SU)>0.1){
#   $7="PASS"; print $0;
#  } else if ( I$SVTYPE=="BND" && $8>=100 && (I$SR/I$SU)>=0.25 && (I$PE/I$SU)>=0.25){
#   $7="PASS"; print $0;
#  } else {
#   $7="LOW"; print $0;
#  }
#}' | bgzip -c > reclassed.filtered.vcf.gz

```
