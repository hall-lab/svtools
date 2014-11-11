SVTools
========
Tools for processing and analyzing structural variants.

## Table of contents
* [bedpeToVcf](#bedpetovcf)

### bedpeToVcf

Convert a LUMPY bedpe file to VCF

#### Usage
```
usage: bedpeToVcf [-h] [-t TOOL] -c SAMPLE_CONFIG [-f FASTA] [-b BEDPE]
                  [-o OUTPUT]

options:
  -h, --help            show this help message and exit
  -t TOOL, --tool TOOL  Tool used to generate calls
  -c SAMPLE_CONFIG, --sample_config SAMPLE_CONFIG
                        Tab delimited sample config file of NAME id TYPE
                        (Example: NA12878 10  PE)
  -f FASTA, --fasta FASTA
                        Indexed fasta file of the reference genome
  -b BEDPE, --bedpe BEDPE
                        BEDPE input (default: stdin)
  -o OUTPUT, --output OUTPUT
                        Output VCF to write (default: stdout)
```


#### Example
```
bedpeToVcf -b samples.sv.bedpe -c samples.config
```
```
##fileformat=VCFv4.1
##fileDate=20141111
##reference=
##INFO=<ID=TOOL,Number=1,Type=String,Description="Tool used to generate variant call">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=STR,Number=.,Type=String,Description="Strand orientation of the adjacency in BEDPE format">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=BKPTID,Number=.,Type=String,Description="ID of the assembled alternate allele in the assembly file">
##INFO=<ID=PARID,Number=1,Type=String,Description="ID of partner breakend">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=SUP,Number=.,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">
##INFO=<ID=PESUP,Number=.,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">
##INFO=<ID=SRSUP,Number=.,Type=Integer,Description="Number of split reads supporting the variant across all samples">
##INFO=<ID=EVTYPE,Number=.,Type=String,Description="Type of LUMPY evidence contributing to the variant call">
##INFO=<ID=PRIN,Number=0,Type=Flag,Description="Indicates variant as the principal variant in a BEDPE pair">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=CNV,Description="Copy number variable region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SUP,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">
##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of split reads supporting the variant">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
##FORMAT=<ID=CNL,Number=.,Type=Float,Description="Copy number genotype likelihood form imprecise events">
##FORMAT=<ID=NQ,Number=1,Type=Integer,Description="Phred style probability score that the variant is novel">
##FORMAT=<ID=HAP,Number=1,Type=Integer,Description="Unique haplotype identifier">
##FORMAT=<ID=AHAP,Number=1,Type=Integer,Description="Unique identifier of ancestral haplotype">
#CHROM							    POS	       ID REF	    ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1							    869423     1  N	    <DEL>	0.00	.	SVTYPE=DEL;SVLEN=-857;END=870280;STR=+-:25;IMPRECISE;CIPOS=-1,34;CIEND=0,0;EVENT=1;SUP=25;PESUP=25;SRSUP=0;EVTYPE=PE;PRIN GT:SUP:PE:SR 0/1:25:25:0
1							    1531294    2  N	    <DUP>	0.00	.	SVTYPE=DUP;SVLEN=344;END=1531638;STR=-+:4;IMPRECISE;CIPOS=-138,1;CIEND=-1,174;EVENT=2;SUP=4;PESUP=4;SRSUP=0;EVTYPE=PE;PRIN	       GT:SUP:PE:SR 0/1:4:4:0
1							    1530938    3  N	    <DEL>	0.00	.	SVTYPE=DEL;SVLEN=-595;END=1531533;STR=+-:4;IMPRECISE;CIPOS=-1,165;CIEND=0,0;EVENT=3;SUP=4;PESUP=4;SRSUP=0;EVTYPE=PE;PRIN	       GT:SUP:PE:SR 0/1:4:4:0
1							    1584545    4  N	    <DEL>	0.00	.	SVTYPE=DEL;SVLEN=-63155;END=1647700;STR=+-:4;IMPRECISE;CIPOS=0,180;CIEND=0,0;EVENT=4;SUP=4;PESUP=4;SRSUP=0;EVTYPE=PE;PRIN	       GT:SUP:PE:SR 0/1:4:4:0
1							    1588585    5  N	    <DUP>	0.00	.	SVTYPE=DUP;SVLEN=65356;END=1653941;STR=-+:7;IMPRECISE;CIPOS=-126,1;CIEND=-2,67;EVENT=5;SUP=7;PESUP=7;SRSUP=0;EVTYPE=PE;PRIN	       GT:SUP:PE:SR 0/1:7:7:0
1							    1594964    6  N	    <DUP>	0.00	.	SVTYPE=DUP;SVLEN=65855;END=1660819;STR=-+:8;IMPRECISE;CIPOS=-81,2;CIEND=-1,127;EVENT=6;SUP=8;PESUP=8;SRSUP=0;EVTYPE=PE;PRIN	       GT:SUP:PE:SR 0/1:8:8:0
1							    2566176    7  N	    <DEL>	0.00	.	SVTYPE=DEL;SVLEN=-418;END=2566594;STR=+-:14;IMPRECISE;CIPOS=-2,68;CIEND=0,0;EVENT=7;SUP=14;PESUP=14;SRSUP=0;EVTYPE=PE;PRIN	       GT:SUP:PE:SR 0/1:14:14:0
1							    2911548    8  N	    <DEL>	0.00	.	SVTYPE=DEL;SVLEN=-302;END=2911850;STR=+-:20;CIPOS=0,0;CIEND=0,0;EVENT=8;SUP=20;PESUP=8;SRSUP=12;EVTYPE=PE,SR;PRIN		       GT:SUP:PE:SR 0/1:20:8:12
1							    2919034    9  N	    <DEL>	0.00	.	SVTYPE=DEL;SVLEN=-332;END=2919366;STR=+-:22;CIPOS=0,0;CIEND=0,0;EVENT=9;SUP=22;PESUP=10;SRSUP=12;EVTYPE=PE,SR;PRIN		       GT:SUP:PE:SR 0/1:22:10:12
1							    3092611    10 N	    <DUP>	0.00	.	SVTYPE=DUP;SVLEN=251;END=3092862;STR=-+:4;IMPRECISE;CIPOS=-3,1;CIEND=-3,3;EVENT=10;SUP=4;PESUP=3;SRSUP=1;EVTYPE=PE,SR;PRIN	       GT:SUP:PE:SR 0/1:4:3:1
1							    4124445    11 N	    <DUP>	0.00	.	SVTYPE=DUP;SVLEN=550;END=4124995;STR=-+:6;IMPRECISE;CIPOS=-1,1;CIEND=0,0;EVENT=11;SUP=6;PESUP=4;SRSUP=2;EVTYPE=PE,SR;PRIN	       GT:SUP:PE:SR 0/1:6:4:2
1							    4124747    12 N	    <DEL>	0.00	.	SVTYPE=DEL;SVLEN=-1589;END=4126336;STR=+-:4;IMPRECISE;CIPOS=-2,3;CIEND=-130,1;EVENT=12;SUP=4;PESUP=4;SRSUP=0;EVTYPE=PE;PRIN	       GT:SUP:PE:SR 0/1:4:4:0
1							    4125505    13 N	    <DEL>	0.00	.	SVTYPE=DEL;SVLEN=-1554;END=4127059;STR=+-:5;IMPRECISE;CIPOS=-1,165;CIEND=0,0;EVENT=13;SUP=5;PESUP=5;SRSUP=0;EVTYPE=PE;PRIN	       GT:SUP:PE:SR 0/1:5:5:0
1							    5447229    14 N	    <DUP>	0.00	.	SVTYPE=DUP;SVLEN=210;END=5447439;STR=-+:11;CIPOS=0,0;CIEND=0,0;EVENT=14;SUP=11;PESUP=1;SRSUP=10;EVTYPE=PE,SR;PRIN		       GT:SUP:PE:SR 0/1:11:1:10
```

Example sample config file (tab delimited)
```
NA12877	10	PE
NA12877	11	SR
NA12877	12	RD
NA12878	20	PE
NA12878	21	SR
NA12878	22	RD
NA12879	30	PE
NA12879	31	SR
NA12879	32	RD
NA12880	40	PE
NA12880	41	SR
NA12880	42	RD
```
