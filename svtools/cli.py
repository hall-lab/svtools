import argparse, sys
import svtools.lsort
import svtools.lmerge
import svtools.vcfpaste
import svtools.copynumber
import svtools.afreq
import svtools.bedpetobed12
import svtools.bedpetovcf
import svtools.vcftobedpe
import svtools.vcfsort
import svtools.bedpesort
import svtools.genotype

def svtools_cli_parser():
    parser = argparse.ArgumentParser(description='Comprehensive utilities to explore structural variation in genomes', prog='svtools')
    subparsers = parser.add_subparsers(title=None, metavar='subcommand', help='description')

    lsort = subparsers.add_parser('lsort', help=svtools.lsort.description())
    svtools.lsort.add_arguments_to_parser(lsort)

    lmerge = subparsers.add_parser('lmerge', help=svtools.lmerge.description())
    svtools.lmerge.add_arguments_to_parser(lmerge)
    
    vcf_paste = subparsers.add_parser('vcfpaste', help=svtools.vcfpaste.description())
    svtools.vcfpaste.add_arguments_to_parser(vcf_paste)

    copynumber = subparsers.add_parser('copynumber', help=svtools.copynumber.description())
    svtools.copynumber.add_arguments_to_parser(copynumber)

    genotype = subparsers.add_parser('genotype', help=svtools.genotype.description())
    svtools.genotype.add_arguments_to_parser(genotype)
    
    afreq = subparsers.add_parser('afreq', help=svtools.afreq.description())
    svtools.afreq.add_arguments_to_parser(afreq)

    bedpetobed12 = subparsers.add_parser('bedpetobed12', help=svtools.bedpetobed12.description())
    svtools.bedpetobed12.add_arguments_to_parser(bedpetobed12)

    bedpetovcf = subparsers.add_parser('bedpetovcf', help=svtools.bedpetovcf.description())
    svtools.bedpetovcf.add_arguments_to_parser(bedpetovcf)

    vcftobedpe = subparsers.add_parser('vcftobedpe', help=svtools.vcftobedpe.description())
    svtools.vcftobedpe.add_arguments_to_parser(vcftobedpe)

    vcfsort = subparsers.add_parser('vcfsort', help=svtools.vcfsort.description())
    svtools.vcfsort.add_arguments_to_parser(vcfsort)

    bedpesort = subparsers.add_parser('bedpesort', help=svtools.bedpesort.description())
    svtools.bedpesort.add_arguments_to_parser(bedpesort)
    
    return parser

if __name__ == '__main__':
    parser = svtools_cli_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))


