#!/usr/bin/bash

# This hits files on our filesystem to create our test files.

zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12878/NA12878.sv.vcf.gz | head -1000 > NA12878.vcf
