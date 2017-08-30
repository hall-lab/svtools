# svtools - Comprehensive utilities to explore structural variations in genomes

[![License](https://img.shields.io/github/license/hall-lab/svtools.svg)](LICENSE.txt)
[![Build Status](https://travis-ci.org/hall-lab/svtools.svg?branch=master)](https://travis-ci.org/hall-lab/svtools) 
[![Coverage Status](https://coveralls.io/repos/github/hall-lab/svtools/badge.svg?branch=master)](https://coveralls.io/github/hall-lab/svtools?branch=master)

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![PyPI](https://img.shields.io/pypi/v/svtools.svg)](https://pypi.python.org/pypi/svtools)
[![DOI](https://zenodo.org/badge/16104/hall-lab/svtools.svg)](https://zenodo.org/badge/latestdoi/16104/hall-lab/svtools)

## Summary
`svtools` is a suite of utilities designed to help bioinformaticians construct and explore cohort-level structural variation calls. It is designed to efficiently merge and genotype calls from [`speedseq sv`](https://github.com/hall-lab/speedseq) across thousands to tens of thousands of genomes.

## Table of Contents
1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Obtaining help](#obtaining-help)
4. [Usage](#usage)
5. [Citing `svtools`](#citing-svtools)
6. [Troubleshooting](#troubleshooting)

## Requirements
* A Linux-like environment with bash, awk, and sort
* A working installation of [`cnvnator-multi`](https://github.com/hall-lab/speedseq#cnvnator) (a component of [`speedseq sv`](https://github.com/hall-lab/speedseq))
* [Python2.7](https://www.python.org/)
   * [Numpy](http://www.numpy.org/)
   * [Scipy](https://www.scipy.org/)
   * [Pandas](http://pandas.pydata.org/)
   * [Statsmodels](http://statsmodels.sourceforge.net/)
   * [Pysam](https://github.com/pysam-developers/pysam) (≥0.8.4)
 
## Installation
We recommend you install using `conda`, but you may also install via `pip`. For more detailed instructions, see our [Installation guide](INSTALL.md).

### Installing via conda
```
conda install -c bioconda svtools
```

### Installing via pip
```
pip install svtools
```

## Obtaining help
Please see the documentation on, or linked to, this page. For additional help or to report a bug, please open an issue in the `svtools` repository: https://github.com/hall-lab/svtools/issues

## Usage
`svtools` consists of subcommands for processing VCF or BEDPE files of structural variants and one accessory script (`create_coordinates`).

```
usage: svtools [-h] [--version] [--support] subcommand ...

Comprehensive utilities to explore structural variation in genomes

optional arguments:
  -h, --help     show this help message and exit
  --version      show program's version number and exit
  --support      information on obtaining support

  subcommand     description
    lsort        sort N LUMPY VCF files into a single file
    lmerge       merge LUMPY calls inside a single file from svtools lsort
    vcfpaste     paste VCFs from multiple samples
    copynumber   add copynumber information using cnvnator-multi
    genotype     compute genotype of structural variants based on breakpoint depth
    afreq        add allele frequency information to a VCF file
    bedpetobed12 convert a BEDPE file to BED12 format for viewing in IGV or the
                 UCSC browser
    bedpetovcf   convert a BEDPE file to VCF
    vcftobedpe   convert a VCF file to a BEDPE file
    vcfsort      sort a VCF file
    bedpesort    sort a BEDPE file
    prune        cluster and prune a BEDPE file by position based on allele
                 frequency
    varlookup    look for variants common between two BEDPE files
    classify     reclassify DEL and DUP based on read depth information
```

## Citing svtools
Until `svtools` is published, please cite using its [DOI](https://zenodo.org/badge/latestdoi/16104/hall-lab/svtools). Note that this link corresponds to the latest version. If you used an earlier version then your DOI may be different and you can find it on [Zenodo](https://zenodo.org/search?ln=en&cc=software&p=svtools&action_search=).

## Troubleshooting
As issues arise and common problems are identified, we will list them here.

Note: For additional information and usage refer to the [Tutorial.md](Tutorial.md) file.
