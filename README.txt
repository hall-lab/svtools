===========
SVTOOls
===========
Warning: This suite is under active development and not a stable version. We can't  provide the support for this version of svtools.

svtools: comprehensive utilities to explore structural variations in genomes.

	Hall-lab
	$Revision: 0.0.1 $
	$Date: 2015-10-14 14:31 $

usage:    svtools <subcommand> [options]
The svtools sub-commands include:
[ general utilities ]
  vcftobedpe      converts vcf file into bedpe.
  bedpetovcf      converts bedpe file to vcf.
  bedpetobed12    converts bedpe file to bed12.
  vcfsort         sorts a vcf file.
  bedpesort       sorts a bedpe file.


[ callset generation ]
  prune           cluster a BEDPE file by position based on allele frequency.
  varlookup       look for variants common between two bedpe files.
  afreq           add allele frequency information to a VCF file.
  lsort           sorts a vcf file by type.
  lmerge          merges multiple sorted vcf files.
  genotype        return a vcf file with genotype information added by svtyper.
  copynumber      add cn information using cnvnator.
  vcfpaste        combine multiple vcf files produced by genotype command.
  classify        classify structural variants


[ General help ]
  --help          print this help menu.
  --version       what version of svtools are you using?.
  --contact       feature requests, bugs, mailing lists, etc.


Url <https://github.com/hall-lab/svtools>`.
Note: For additional information and usage refer "svtools_demo.sh" file
