import argparse, sys, errno
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
import svtools.prune
import svtools.varlookup
import svtools.sv_classifier

class SupportAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        support_string = 'For further help or to report a bug, please open an issue on the svtools repository: https://github.com/hall-lab/svtools/issues'
        print support_string
        sys.exit()

def svtools_cli_parser():
    parser = argparse.ArgumentParser(description='Comprehensive utilities to explore structural variation in genomes', prog='svtools')
    version_string = '%(prog)s {0}'.format(svtools.__version__)
    parser.add_argument('--version', action='version', version=version_string)
    parser.add_argument('--support', action=SupportAction, nargs=0, help='information on obtaining support')
    subparsers = parser.add_subparsers(title=None, metavar='subcommand', help='description')

    lsort = subparsers.add_parser('lsort', help=svtools.lsort.description(), epilog=svtools.lsort.epilog())
    svtools.lsort.add_arguments_to_parser(lsort)

    lmerge = subparsers.add_parser('lmerge', help=svtools.lmerge.description(), epilog=svtools.lmerge.epilog())
    svtools.lmerge.add_arguments_to_parser(lmerge)

    vcf_paste = subparsers.add_parser('vcfpaste', help=svtools.vcfpaste.description(), epilog=svtools.vcfpaste.epilog())
    svtools.vcfpaste.add_arguments_to_parser(vcf_paste)

    copynumber = subparsers.add_parser('copynumber', help=svtools.copynumber.description(), epilog=svtools.copynumber.epilog())
    svtools.copynumber.add_arguments_to_parser(copynumber)

    genotype = subparsers.add_parser('genotype', help=svtools.genotype.description())
    svtools.genotype.add_arguments_to_parser(genotype)

    afreq = subparsers.add_parser('afreq', help=svtools.afreq.description(), epilog=svtools.afreq.epilog())
    svtools.afreq.add_arguments_to_parser(afreq)

    bedpetobed12 = subparsers.add_parser('bedpetobed12', help=svtools.bedpetobed12.description(), epilog=svtools.bedpetobed12.epilog())
    svtools.bedpetobed12.add_arguments_to_parser(bedpetobed12)

    bedpetovcf = subparsers.add_parser('bedpetovcf', help=svtools.bedpetovcf.description(), epilog=svtools.bedpetovcf.epilog())
    svtools.bedpetovcf.add_arguments_to_parser(bedpetovcf)

    vcftobedpe = subparsers.add_parser('vcftobedpe', help=svtools.vcftobedpe.description(), epilog=svtools.vcftobedpe.epilog())
    svtools.vcftobedpe.add_arguments_to_parser(vcftobedpe)

    vcfsort = subparsers.add_parser('vcfsort', help=svtools.vcfsort.description())
    svtools.vcfsort.add_arguments_to_parser(vcfsort)

    bedpesort = subparsers.add_parser('bedpesort', help=svtools.bedpesort.description())
    svtools.bedpesort.add_arguments_to_parser(bedpesort)

    prune = subparsers.add_parser('prune', help=svtools.prune.description(), epilog=svtools.prune.epilog())
    svtools.prune.add_arguments_to_parser(prune)

    varlookup = subparsers.add_parser('varlookup', help=svtools.varlookup.description())
    svtools.varlookup.add_arguments_to_parser(varlookup)

    classifier = subparsers.add_parser('classify', help=svtools.sv_classifier.description())
    svtools.sv_classifier.add_arguments_to_parser(classifier)

    return parser

def main():
    parser = svtools_cli_parser()
    args = parser.parse_args()
    try:
        sys.exit(args.entry_point(args))
    except IOError as e:
        if e.errno == errno.EPIPE:
            sys.exit(141)

if __name__ == '__main__':
    main()

