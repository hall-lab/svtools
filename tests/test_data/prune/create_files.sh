#!/usr/bin/bash
set -o pipefail

# This hits files on our filesystem to create our test files.

zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12878/NA12878.sv.vcf.gz | head -1000 > NA12878.vcf
zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12891/NA12891.sv.vcf.gz | head -1000 > NA12891.vcf
zcat /gscmnt/gc2802/halllab/sv_aggregate/MISC/lumpy/NA12892/NA12892.sv.vcf.gz | head -1000 > NA12892.vcf

DIR=`pwd`
echo -e "$DIR/NA12878.vcf\n$DIR/NA12891.vcf\n$DIR/NA12892.vcf" > vcf_list

export PYTHONPATH=~/src/svtools/:$PYTHONPATH 
python ~/src/svtools/svtools/vcfpaste.py -f vcf_list | python ~/src/svtools/svtools/afreq.py | python ~/src/svtools/svtools/vcftobedpe.py > input.bed
cat input.bed | grep -v 'MISSING' > input.no_missing.bed 

rm -f *vcf*
