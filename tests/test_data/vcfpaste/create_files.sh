#!/usr/bin/bash

# This hits files on our filesystem to create our test files.

zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12878/NA12878.sv.vcf.gz | head -50 > NA12878.vcf
bgzip -c NA12878.vcf > NA12878.vcf.gz
zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12891/NA12891.sv.vcf.gz | head -50 > NA12891.vcf
bgzip -c NA12891.vcf > NA12891.vcf.gz
zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12892/NA12892.sv.vcf.gz | head -50 > NA12892.vcf
bgzip -c NA12892.vcf > NA12892.vcf.gz
