#!/usr/bin/env bash

# # make truth set and sort the BAMs
# BAM=/gscmnt/gc2802/halllab/sv_aggregate/MISC/realigned_BAMs/NA12878/NA12878.bam
# BAM_BASE=`basename $BAM .bam`
# ../svtyper \
#     -i example.vcf \
#     -B $BAM \
#     -l $BAM_BASE.bam.json \
#     -w $BAM_BASE.target_loci.bam \
#     > example.gt.vcf
#
# sambamba sort NA12878.target_loci.bam

# run test
../svtyper \
    -i example.vcf \
    -B NA12878.target_loci.sorted.bam \
    -l NA12878.bam.json \
    > out.vcf

diff -I '^##fileDate=' example.gt.vcf out.vcf
