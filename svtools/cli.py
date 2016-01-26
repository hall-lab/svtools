import argparse, sys
import svtools.lsort
import svtools.lmerge
import svtools.vcfpaste
import svtools.afreq

def svtools_cli_parser():
    parser = argparse.ArgumentParser(description='Comprehensive utilities to explore structural variation in genomes', prog='svtools')
    subparsers = parser.add_subparsers(title=None, metavar='subcommand', help='description')

    lsort = subparsers.add_parser('lsort', help=svtools.lsort.description())
    svtools.lsort.add_arguments_to_parser(lsort)

    lmerge = subparsers.add_parser('lmerge', help=svtools.lmerge.description())
    svtools.lmerge.add_arguments_to_parser(lmerge)
    
    vcf_paste = subparsers.add_parser('vcfpaste', help=svtools.vcfpaste.description())
    svtools.vcfpaste.add_arguments_to_parser(vcf_paste)
    
    afreq = subparsers.add_parser('afreq', help=svtools.afreq.description())
    svtools.afreq.add_arguments_to_parser(afreq)

    return parser

if __name__ == '__main__':
    parser = svtools_cli_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))


