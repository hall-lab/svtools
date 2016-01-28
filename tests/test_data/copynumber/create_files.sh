#!/usr/bin/bash

# This hits files on our filesystem to create our test files.

zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12878/NA12878.sv.vcf.gz | grep -v 'SVTYPE=BND' | head -50 > NA12878.vcf
zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12891/NA12891.sv.vcf.gz | grep -v 'SVTYPE=BND' | head -50 > NA12891.vcf
echo "NA12878.vcf" >> vcflist
echo "NA12891.vcf" >> vcflist
PYTHONPATH=~/src/svtools/:$PYTHONPATH python ~/src/svtools/svtools/cli.py vcfpaste -f vcflist > input.vcf
rm -f NA12878.vcf NA12891.vcf vcflist
