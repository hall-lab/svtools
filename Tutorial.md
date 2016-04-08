#Example analysis using `svtools`
This tutorial will help you begin to explore the use of `svtools` to analyze an SV VCF.  It will help you to satisfy the computing environment requirements, gather the required genomic data, and walk through basic analysis using `svtools`.
This tutorial includes example commands that you can alter to refer to your sample names.

##Table of Contents
1. Satisfy computing environment requirements
2. Gather genomic data and generate needed helper files
3. Use `svtools` to create a callset
    1. Use vawk to remove homozygous reference variants from SpeedSeq SV VCFs
    2. Use `svtools lsort` to combine and sort variants from multiple samples
    3. Use `svtools lmerge` to merge variants deemed to be identical in the sorted VCF
    4. (Optional) Remove EBV (Epstein-Barr Virus) variants
    5. Use `svtools genotype` to force genotypes for variant positions discovered in other samples
    6. Use `svtools copynumber` to create per-sample copy number annotations based on CNVnator histograms 
        1. Prepare environment for CNVnator
        2. Make an uncompressed copy 
        3. Make coordinate file
        4. Annotate variants with copy number from CNVnator using `svtools copynumber`
    7. Use `svtools vcfpaste` to construct a VCF that pastes in genotype and copy number information
    8. Use `svtools prune` to filter out additional variants deemed to be identical  

## 1) Satisfy computing environment requirements
### 1) Install SpeedSeq and dependencies
Installation instructions have been provided in the [SpeedSeq github repository](https://github.com/hall-lab/speedseq).
### 2) Install `svtools`
Installation instructions have been provided in the [INSTALL.md](https://github.com/jeldred/svtools/blob/install_documentation/INSTALL.md) and DEVELOPER.md of this repo.

## 2) Gather genomic data and generate needed helper files
### 1) Get or Create SpeedSeq aligned BAMs, splitter files, and SV VCF files
To get a small set of bam files suitable for this tutorial. I recommend getting 3 bams from http://www.ebi.ac.uk/ena/data/view/ERP001960 in the NA12878 pedigree. A simple command line to achieve this is listed below. Or you could try the [Aspera client.](http://downloads.asperasoft.com/connect2/)
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/NA12877_S1.bam
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/NA12878_S1.bam
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/NA12879_S1.bam
```
Downloading these bams will consume a significant amount of time, bandwidth and disk space.
These three files sum to 317GB.
Use documentation on the SpeedSeq github page to produce the required files for the rest of this tutorial at [SpeedSeq github repository](https://github.com/hall-lab/speedseq).
### 2) Get Reference FASTA
We recommend using the GRCh37 human genome for SpeedSeq, available here:
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai
The genome FASTA file should be unzipped and indexed with BWA before running SpeedSeq.
### 3) Create cn.list file
The cn.list file has a single column that contains the path to the VCF files output in the Copy Number Annotation step of this tutorial.

## 3) Use `svtools` to create a callset
### 1) Use vawk to remove homozygous reference variants from SpeedSeq SV VCFs
This step will remove variants that have been detected by Lumpy but then determined to be homozygous reference when SVTyper is run.
This command will need to be run once per sample and ouputs one non_ref VCF file per sample.
```
  zcat NA12877.sv.vcf.gz \
  | vawk --header '{if(S$*$GT!="0/0" && S$*$GT!="./.") print $0}' \
  > NA12877.sv.non_ref.vcf
```

### 2) Use `svtools lsort` to combine and sort variants from multiple samples
`svtools lsort` takes a space separated list of all of your non_ref vcf files as arguments.
The example below shows us combining three samples.  The output of this step is one sorted and compressed VCF file containing all variants detected in the three input files.
```
svtools lsort NA12877.sv.non_ref.vcf NA12878.sv.non_ref.vcf NA12879.sv.non_ref.vcf \
| bgzip -c > sorted.vcf.gz
```
```
:note svtools lsort will remove variants with the SECONDARY tag in the INFO field.
This will cause the sorted VCF to have fewer variant lines than the input.
```
### 3) Use `svtools lmerge` to merge variants deemed to be identical in the sorted VCF
```
zcat sorted.vcf.gz \
| svtools lmerge -i /dev/stdin --product -f 20 \
| bgzip -c > merged.vcf.gz "

```
```
:note svtools lmerge will return variant lines for SECONDARY break ends in addition to merging variants.
This will sometimes cause the merged VCF to have more variant lines than the input.
```
### 4) (Optional) Remove variants detected by alignment to the EBV (Epstein-Barr Virus) contig
If your reference contains a contig representing EBV then you may wish to remove SVs involved with this sequence.
```
zcat merged.vcf.gz \
| vawk --header '{if($0 !~ /NC_007605/) print $0}' \
| bgzip -c > merged.no_EBV.vcf.gz
```

### 5) Use `svtools genotype` to force genotypes for variant positions discovered in other samples
`svtools genotype` will calculate a genotype for each sample at the variant positions in the merged.no_EBV.vcf.gz file.
It requires the aligned BAM and a splitters BAM file for each sample. This step will output a fully genotyped VCF file for each sample.
You will also need to prepare a gt subdirectory to store the output of these commands to avoid name collisions with the upcoming copy number output.
```
mkdir -p gt

"zcat merged.no_EBV.vcf.gz \
| vawk --header '{  \$6=\".\"; print }' \
| svtools genotype \
  -B NA12877.bam \
  -S NA12877.splitters.bam \
| sed 's/PR...=[0-9\.e,-]*\(;\)\{0,1\}\(\t\)\{0,1\}/\2/g' - \
> gt/NA12877.vcf"
```
You will need to repeat the `svtools genotype` command above for the other two samples (NA12878,NA12879) as well.

### 6) Use `svtools copynumber` to create per-sample copy number annotations based on CNVnator histograms 
#### 1) Prepare environment for CNVnator
CNVnator require the ROOT package to function. This file must be sourced before running CNVnator.
```
source /gsc/pkg/root/root/bin/thisroot.sh
```
#### 2) Make an uncompressed copy 
```
zcat merged.no_EBV.vcf.gz > merged.no_EBV.vcf
```
#### 3) Make coordinate file
CNVnator will return the copynumber for a list of coordinates. This script will create such a list and is deployed upon installation of `svtools`.
```
create_coordinates -i merged.no_EBV.vcf -o coordinates
```

#### 4) Annotate variants with copy number from CNVnator using `svtools copynumber`
This step assumes you have already run CNVnator and that the output required for this step is stored in your analysis directory at 
`/temp/cnvnator-temp/NA12877.bam.hist.root`. If you have installed SpeedSeq, CNVnator is run as part of `speedseq sv`. More details about `speedseq sv` are [here](https://github.com/hall-lab/speedseq#speedseq-sv)

You will also need to prepare a subdirectory to hold the Copy Number(cn) vcf files 
```
mkdir -p cn
```

Then run `svtools copynumber` to add in copynumber values to non-BND variants.
```
svtools copynumber \
--cnvnator cnvnator-multi \
-s NA12877 \
-w 100 \
-r /temp/cnvnator-temp/NA12877.bam.hist.root \
 -c coordinates \
 -v gt/NA12877.vcf \
> cn/NA12877.vcf
```

### 7) Use `svtools vcfpaste` to construct a VCF that pastes in genotype and copy number information
`svtools vcfpaste` takes the list of the VCFs generated that contain the additional information for every sample that we have been building up step by step.  In this tutorial we call that file cn.list and it contains one column that holds the path to the VCF files generated in the previous step.
```
svtools vcfpaste \
-m merged.no_EBV.vcf \
-f cn.list \
-q \
| bgzip -c \
> merged.sv.gt.cn.vcf.gz
```

### 8) Use `svtools prune` to filter out additional variants deemed to be identical
```
bsub -q long -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' "zcat merged.sv.gt.cn.vcf.gz \
| svtools afreq \
| svtools vcftobedpe \
| svtools bedpesort \
| svtools prune -s -d 100 -e \"AF\" \
| svtools bedpetovcf \
| bgzip -c > merged.sv.new_pruned.vcf.gz"
```

