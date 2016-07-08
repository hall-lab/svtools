#!/usr/bin/bash
set -o pipefail

# This hits files on our filesystem to create our test files.

zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12878/NA12878.sv.vcf.gz | grep -v '0\/0' > NA12878.vcf
zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12891/NA12891.sv.vcf.gz | grep -v '0\/0' > NA12891.vcf
zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12892/NA12892.sv.vcf.gz | grep -v '0\/0' > NA12892.vcf

#svtools lsort *vcf > lsort_expected
#rm *vcf
