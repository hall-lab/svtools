import argparse, sys
import svtools.lsort
import svtools.vcfpaste
import svtools.afreq

def svtools_cli_parser():
    parser = argparse.ArgumentParser(description='Comprehensive utilities to explore structural variation in genomes', prog='svtools')
    subparsers = parser.add_subparsers(title='subcommands',
                                       description='valid subcommands',
                                       help='additional help')

    lsort = subparsers.add_parser('lsort', description=svtools.lsort.description())
    svtools.lsort.add_arguments_to_parser(lsort)
    
    vcf_paste = subparsers.add_parser('vcfpaste', description=svtools.vcfpaste.description())
    svtools.vcfpaste.add_arguments_to_parser(vcf_paste)
    
    afreq = subparsers.add_parser('afreq', description=svtools.afreq.description())
    svtools.afreq.add_arguments_to_parser(afreq)

    return parser

if __name__ == '__main__':
    parser = svtools_cli_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))


