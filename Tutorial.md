#Example analysis using `svtools`
This tutorial will help you begin to explore the use of `svtools` to analyze an SV VCF.  It will help you to satisfy the computing environment requirements, gather the required genomic data, and walk through basic analysis using `svtools`.
This tutorial includes example commands that you can alter to refer to your sample names.

##Table of Contents
1. Satisfy computing environment requirements
2. Gather genomic data and generate needed helper files
3. Use `svtools` to create a callset
    1. Use vawk to remove homozygous reference variants from SpeedSeq SV VCFs
    2. Use `svtools lsort` to combine and sort variants from multiple samples
    3. Use `svtools lmerge` to merge variant calls likely representing the same variant in the sorted VCF
    4. Use `svtools genotype` to genotype all samples for all variants present in the merged set for variant positions discovered in other samples
    5. Use `svtools copynumber` to create per-sample copynumber annotations based on CNVnator histograms 
        1. Prepare environment for CNVnator
        2. Make an uncompressed copy 
        3. Make coordinate file
        4. Annotate variants with copynumber from CNVnator using `svtools copynumber`
    6. Use `svtools vcfpaste` to construct a VCF that pastes together the individual genotyped and copynumber annotated vcfs
    7. Use `svtools prune` to filter out additional variant calls likely representing the same variant  

## Satisfy computing environment requirements
### Install SpeedSeq and dependencies
Installation instructions for SpeedSeq are available in the [SpeedSeq github repository](https://github.com/hall-lab/speedseq).
### Install `svtools`
Installation instructions have been provided in the [INSTALL.md](INSTALL.md) and [DEVELOPER.md](DEVELOPER.md) of this repo.

## Gather genomic data and generate needed helper files

### Get a Reference FASTA
We recommend using the GRCh37 human genome for SpeedSeq, available here:
* ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
* ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai

**Note:** The genome FASTA file should be unzipped and indexed with BWA before running SpeedSeq.

### Get or Create SpeedSeq aligned BAMs, splitter BAM files, and SV VCF files
To get a small set of BAM files suitable for this tutorial, we recommend getting 3 BAMs from http://www.ebi.ac.uk/ena/data/view/ERP001960 in the NA12878 pedigree. A simple command line to achieve this is listed below. Or you could try the [Aspera client.](http://downloads.asperasoft.com/connect2/)
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/NA12877_S1.bam
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/NA12878_S1.bam
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/NA12879_S1.bam
```
Downloading these BAMs will consume a significant amount of time, bandwidth and disk space (~317GB).
Follow the documentation on the [SpeedSeq Github page](https://github.com/hall-lab/speedseq) to run `speedseq realign` and `speedseq sv` on these BAMs. This will produce the required files for the rest of this tutorial.

## Use `svtools` to create a callset
### Use vawk to remove homozygous reference variants from SpeedSeq SV VCFs
This step will remove variants that have been detected by Lumpy but then determined to be homozygous reference when SVTyper is run.
This command will need to be run once per sample and ouputs one non_ref VCF file per sample.
```
  zcat NA12877.sv.vcf.gz \
  | vawk --header '{if(S$*$GT!="0/0" && S$*$GT!="./.") print $0}' \
  > NA12877.sv.non_ref.vcf
```
**Note:** vawk is included as part of SpeedSeq.

### Use `svtools lsort` to combine and sort variants from multiple samples
`svtools lsort` takes a space separated list of all of the non_ref VCF files generated in the previous step as arguments.
The example below shows us combining three samples.  The output of this step is one sorted and compressed VCF file containing all variants detected in the three input files.
```
svtools lsort NA12877.sv.non_ref.vcf NA12878.sv.non_ref.vcf NA12879.sv.non_ref.vcf \
| bgzip -c > sorted.vcf.gz
```

**Note:** `svtools lsort` will remove variants with the SECONDARY tag in the INFO field.
This will cause the sorted VCF to have fewer variant lines than the input.

###Use `svtools lmerge` to merge variant calls likely representing the same variant in the sorted VCF
```
zcat sorted.vcf.gz \
| svtools lmerge -i /dev/stdin --product -f 20 \
| bgzip -c > merged.vcf.gz "

```

**Note:** `svtools lmerge` will return variant lines for SECONDARY breakends in addition to merging variants.
This will sometimes cause the merged VCF to have more variant lines than the input.

### Use `svtools genotype` to genotype all samples for all variants present in the merged set
`svtools genotype` will calculate a genotype for each sample at the variant positions in the merged.vcf.gz file.
It requires the aligned BAM and a splitters BAM file for each sample. This step will output a fully genotyped VCF file for each sample.
You will also need to prepare a gt subdirectory to store the output of these commands to avoid name collisions with the upcoming copynumber output.
```
mkdir -p gt

"zcat merged.vcf.gz \
| vawk --header '{  \$6=\".\"; print }' \
| svtools genotype \
  -B NA12877.bam \
  -S NA12877.splitters.bam \
| sed 's/PR...=[0-9\.e,-]*\(;\)\{0,1\}\(\t\)\{0,1\}/\2/g' - \
> gt/NA12877.vcf"
```
You will need to repeat the `svtools genotype` command above for the other two samples (NA12878, NA12879) as well.

### Use `svtools copynumber` to create per-sample copynumber annotations based on CNVnator histograms 
#### Prepare environment for CNVnator
CNVnator requires the ROOT package to function. This file must be sourced before running CNVnator.
```
source /gsc/pkg/root/root/bin/thisroot.sh
```
#### Make an uncompressed copy
This will be used to create the coordinate file in the next step.
```
zcat merged.vcf.gz > merged.vcf
```
#### Make coordinate file
CNVnator will return the copynumber for a list of coordinates. This script will create such a list and is deployed upon installation of `svtools`.
```
create_coordinates -i merged.vcf -o coordinates
```
**Note:** The last line of this file should be the word "exit". This is intentional and [required by CNVnator](https://github.com/abyzovlab/CNVnator).
.
#### Annotate variants with copynumber from CNVnator using `svtools copynumber`
This step assumes you have already run CNVnator and that the output required for this step is stored in your analysis directory at 
`/temp/cnvnator-temp/NA12877.bam.hist.root`. If you have installed SpeedSeq, CNVnator is run as part of `speedseq sv`. More details about `speedseq sv` are [here](https://github.com/hall-lab/speedseq#speedseq-sv)

You will also need to prepare a subdirectory to hold the copynumber(cn) VCF files 
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
**Note:** The argument to the `--cnvnator` option of `svtools copynumber` may need to be the full path to the cnvnator-multi executable included as part of SpeedSeq. This example assumes cnvnator-multi is installed system-wide. 

### Use `svtools vcfpaste` to construct a VCF that pastes together the individual genotyped and copynumber annotated vcfs
`svtools vcfpaste` takes the list of the VCFs generated that contain the additional information for every sample that we have been building up step by step.  In this tutorial we call that file cn.list and it contains one column that holds the path to the VCF files generated in the previous step.

To generate the cn.list file:
```
ls -1 cn/*vcf > cn.list
```

Then run `svtools vcfpaste` to re-assemble a cohort-level VCF file
```
svtools vcfpaste \
-m merged.vcf \
-f cn.list \
-q \
| bgzip -c \
> merged.sv.gt.cn.vcf.gz
```

### Use `svtools prune` to filter out additional variants deemed to be identical
```
bsub -q long -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' "zcat merged.sv.gt.cn.vcf.gz \
| svtools afreq \
| svtools vcftobedpe \
| svtools bedpesort \
| svtools prune -s -d 100 -e \"AF\" \
| svtools bedpetovcf \
| bgzip -c > merged.sv.new_pruned.vcf.gz"
```

