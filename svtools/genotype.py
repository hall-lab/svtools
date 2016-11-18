import sys
import argparse
from svtools.external_cmd import ExternalCmd

class GenotypeVariants(ExternalCmd):
    def __init__(self):
        super(GenotypeVariants, self).__init__('svtyper', 'bin/svtyper/svtyper')

    @staticmethod
    def svtyper_option_lut():
        opts = { 
                'input_vcf' : '-i',
                'output_vcf' : '-o',
                'bam' : '-B',
                'lib_info' : '-l',
                'min_aligned' : '-m',
                'num_samp' : '-n',
                'split_weight' : '--split_weight',
                'disc_weight' : '--disc_weight',
                'write_alignment' : '-w'
                }
        return opts

def description():
    return 'compute genotype of structural variants based on breakpoint depth'

def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--input_vcf', help='VCF input (default: stdin)')
    parser.add_argument('-o', '--output_vcf', help='output VCF to write (default: stdout)')
    parser.add_argument('-B', '--bam', type=str, required=True, help='BAM or CRAM file(s), comma-separated if genotyping multiple BAMs')
    parser.add_argument('-l', '--lib_info', type=str, required=False, help='create/read JSON file of library information')
    parser.add_argument('-m', '--min_aligned', type=int, required=False, default=20, help='minimum number of aligned bases to consider read as evidence [20]')
    parser.add_argument('-n', dest='num_samp', type=int, required=False, default=1000000, help='number of pairs to sample from BAM file for building insert size distribution [1000000]')
    parser.add_argument('--split_weight', type=float, required=False, default=1, help='weight for split reads [1]')
    parser.add_argument('--disc_weight', type=float, required=False, default=1, help='weight for discordant paired-end reads [1]')
    parser.add_argument('-w', '--write_alignment', type=str, required=False, default=None, help='write relevant reads to BAM file')
    parser.add_argument('--debug', action='store_true', help=argparse.SUPPRESS)
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    argdict = vars(args)

    genotyper = GenotypeVariants()
    opts = list()
    optlut = genotyper.svtyper_option_lut()
    for variable, value in argdict.iteritems():
        if variable != 'entry_point' and value not in (False, None):
            opts.extend([optlut[variable], str(value)])
    genotyper.run_cmd_with_options(opts)

if __name__ == "__main__":
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
