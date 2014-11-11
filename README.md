SVTools
========
Tools for processing and analyzing structural variants.

## Table of contents
* [bedpeToVcf](#bedpetovcf)
* [vcfToBedpe](#vcftobedpe)
* [lumpyToBedpe](#lumpyToBedpe)

### bedpeToVcf

Convert a LUMPY BEDPE file to VCF

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
                          (pysam required if using this option)
  -b BEDPE, --bedpe BEDPE
                        BEDPE input (default: stdin)
  -o OUTPUT, --output OUTPUT
                        Output VCF to write (default: stdout)
```


#### Example
```
bedpeToVcf -b samples.sv.bedpe -c samples.config > samples.sv.vcf
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

### vcfToBedpe

Convert a VCF file to a BEDPE file

#### Usage
```
usage: vcfToBedpe [-h] [-v VCF] [-o OUTPUT]

vcfToBedpe
author: Colby Chiang (cc2qe@virginia.edu)
version: $Revision: 0.0.1 $
description: Convert a VCF file to a BEDPE file

optional arguments:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     VCF input (default: stdin)
  -o OUTPUT, --output OUTPUT
                        Output BEDPE to write (default: stdout)
```

#### Example
```
vcfToBedpe -v samples.sv.vcf > samples.sv.bedpe
```

### lumpyToBedpe

Convert a LUMPY BEDPE file to an extended BEDPE for
easier filtering

#### Usage
```
usage: lumpyToBedpe [options]

options:
  -h, --help            show this help message and exit
  -b BEDPE_FILE, --bedpe_file=BEDPE_FILE
                        BEDPE file
  -c CONFIG_FILE, --config_file=CONFIG_FILE
                        Tab-delim sample config file of NAME id TYPE.
                        Example:NA12878 10  PE
```

#### Example
```
lumpyToBedpe -b NA12878_S1.sv.bedpe -c NA12878_S1.sample.config > NA12878_S1.sv.ext.bedpe
```
