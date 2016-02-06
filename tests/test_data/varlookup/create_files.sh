#!/usr/bin/bash
set -o pipefail

# This hits files on our filesystem to create our test files.

zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12878/NA12878.sv.vcf.gz | head -1000 > NA12878.vcf
zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12891/NA12891.sv.vcf.gz | head -1000 > NA12891.vcf

export PYTHONPATH=~/src/svtools/:$PYTHONPATH 
python ~/src/svtools/svtools/afreq.py NA12878.vcf | python ~/src/svtools/svtools/vcftobedpe.py | grep -v 'MISSING' > input_a.bed
python ~/src/svtools/svtools/afreq.py NA12891.vcf | python ~/src/svtools/svtools/vcftobedpe.py | grep -v 'MISSING' > input_b.bed

rm -f *vcf*
