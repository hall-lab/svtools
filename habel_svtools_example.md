#Example analysis using `svtools`
This tutorial will help you begin to explore the use of `svtools` to analyze an SV vcf.  It will help you to satisfy the computing environment requirements, gather the required genomic data, and try basic analysis using `svtools`.
This tutorial includes example commands that you can alter to refer to your sample names.

1. Satisfy computing environment requirements
2. Gather genomic data and generate needed helper files
3. Use `svtools` to create a callset
4. Other Tutorial Resources

##Satisfy computing environment requirements
1. Install SpeedSeq and dependencies 
2. Install `svtools`

###Install SpeedSeq and dependencies
[Installation instructions](https://github.com/hall-lab/speedseq/blob/master) have been provided in SpeedSeq github repository.
###Install `svtools`
[Installation instructions](https://github.com/jeldred/svtools/blob/install_documentation/INSTALL.md) have been provided in the INSTALL.md and DEVELOPER.MD of this repo.

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

##Use `svtools` to create a callset
1. Use vawk to remove 'REF' variants from SpeedSeq SV vcf
2. Use `svtools lsort` to combine and sort variants from multiple samples
3. Use `svtools lmerge` to merge variants deemed to be identical in the sorted vcf
4. Remove EBV (Epstein-Barr Virus) variants
5. Use `svtools genotype` to force genotypes for variant positions discovered in other samples.
6. Use `svtools copynumber` to create per sample copy number annotations based on CNVnator histograms 
    1. prepare environemnt for cnvnator
    2. make uncompressed copy 
    3. make coordinate file
    4. Annotate variants with copy number from cnvnator using `svtools copynumber`
7. Use `svtools vcfpaste` to construct a bam that pastes in genotype and copy number information
8. Use `svtools prune` to filter out additional variants deemed to be identical  
9. Training (exclude?)
10. Classification (exclude?)

###Use vawk to remove 'REF' variants from SpeedSeq SV vcf
This step will remove variants that have been detected but then determined to be homozygous reference.
This command will need to be run once per sample and ouput one non_ref vcf file per sample.
```
  zcat SAMPLE1.sv.vcf.gz \
  | vawk --header '{if(S$*$GT!="0/0" && S$*$GT!="./.") print $0}' \
  > SAMPLE1.sv.non_ref.vcf
```

### Use `svtools lsort` to combine and sort variants from multiple samples
`svtools lsort` takes the space sorted list of all of your non_ref vcf files as arguments.
The example shows us combining three samples.  The output of this step is one sorted and compressed vcf file.
```
svtools lsort SAMPLE1.sv.non_ref.vcf SAMPLE2.sv.non_ref.vcf SAMPLE3.sv.non_ref.vcf \
| bgzip -c > sorted.vcf.gz
```
```
:note svtools lsort will remove variants with the SECONDARY tag in the INFO field.
This will cause the sorted vcf to have fewer variant lines than the input.
```
### Use `svtools lmerge` to merge variants deemed to be identical in the sorted vcf
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
| vawk --header '{if($0 !~ /NC_007605/) print $0}' \
| bgzip -c > merged.no_EBV.vcf.gz
```

### Use `svtools genotype` to force genotypes for variant positions discovered in other samples
`svtools genotype` will calculate a genotype for each sample at the variant positions in the merged.no_EBV.vcf.gz file.
It requires the aligned bam and splitters file for each sample. This step will output a fully genotyped vcf file for each sample.
You will also need to prepare a gt subdirectory to store the output of these commands to avoid name collisions with the upcoming copy number output.
```
mkdir -p gt

"zcat merged.no_EBV.vcf.gz \
| vawk --header '{  \$6=\".\"; print }' \
| svtools genotype \
  -B SAMPLE1.bam \
  -S SAMPLE1.splitters.bam \
| sed 's/PR...=[0-9\.e,-]*\(;\)\{0,1\}\(\t\)\{0,1\}/\2/g' - \
> gt/SAMPLE1.vcf"
```

###Use `svtools copynumber` to create per sample copy number annotations based on CNVnator histograms 
#### prepare environemnt for cnvnator
```
source /gsc/pkg/root/root/bin/thisroot.sh
```
#### make uncompressed copy 
```
zcat merged.no_EBV.vcf.gz > merged.no_EBV.vcf
```
#### make coordinate file
```
create_coordinates -i merged.no_EBV.vcf -o coordinates
```
This step assumes you have already run CNVnator and that the output required for this step is stored in your analysis directory at 
```
/temp/cnvnator-temp/SAMPLE1.bam.hist.root
```
more details about running CNVnator are [not actually here](http://some_link.gsc.wustl.edu)

You will also need to prepare a subdirectory to hold the Copy Number(cn) vcf files 
```
mkdir -p cn

"svtools copynumber \
--cnvnator cnvnator-multi \
-s SAMPLE1 \
-w 100 \
-r /temp/cnvnator-temp/SAMPLE1.bam.hist.root \
 -c coordinates \
 -v gt/SAMPLE1.vcf \
> cn/SAMPLE.vcf"
```

### Use `svtools vcfpaste` to construct a bam that pastes in genotype and copy number information
`svtools vcfpaste` takes the list of the cn vcfs that contain the information that we now want to paste to the end of variant lines that we have been building up step by step.  In this tutorial we call that file cn.list and it contains one column that holds the path to those cn vcf files.
```
"svtools vcfpaste \
-m merged.no_EBV.vcf \
-f cn.list \
-q \
| bgzip -c \
> merged.sv.gt.cn.vcf.gz"
```

### Use `svtools prune` to filter out additional variants deemed to be identical
```
:note Future improvements to svtools lmerge may make svtools prune unnecessary
```
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
 | vawk --header '{if(I$FINMETSEQ_HQ_AF>0) print $0}' \
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

##Other Tutorial Resources
###Example Bash Script
Creating a sample.map file and running lumpy_pipeline.sh is another way to investigate the usage 
of svtools for creating a callset.  It has some additional requirements that are documented by comments in the script.  

[lumpy_pipeline.sh](https://github.com/jeldred/svtools/edit/install_documentation/lumpy_pipeline.sh)
