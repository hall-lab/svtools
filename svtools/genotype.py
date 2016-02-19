import sys
import argparse
from svtools.external_cmd import ExternalCmd

class GenotypeVariants(ExternalCmd):
    def __init__(self):
        super(GenotypeVariants, self).__init__('svtyper', 'bin/svtyper/svtyper')

    @staticmethod
    def svtyper_option_lut():
        opts = { 
                'bam' : '-B',
                'split_bam' : '-S',
                'input_vcf' : '-i',
                'output_vcf' : '-o',
                'splflank' : '-f',
                'discflank' : '-F',
                'split_weight' : '--split_weight',
                'disc_weight' : '--disc_weight',
                'num_samp' : '-n',
                'legacy' : '-M',
                'debug' : '--debug',
                }
        return opts

def description():
    return 'compute genotype of structural variants based on breakpoint depth'

def add_arguments_to_parser(parser):
    parser.add_argument('-B', '--bam', type=str, required=True, help='BAM file(s), comma-separated if genotyping multiple BAMs')
    parser.add_argument('-S', '--split_bam', type=str, required=False, help='split-read bam file for sample, comma-separated if genotyping multiple BAMs')
    parser.add_argument('-i', '--input_vcf', help='VCF input (default: stdin)')
    parser.add_argument('-o', '--output_vcf', help='output VCF to write (default: stdout)')
    parser.add_argument('-f', '--splflank', type=int, required=False, default=20, help='min number of split read query bases flanking breakpoint on either side [20]')
    parser.add_argument('-F', '--discflank', type=int, required=False, default=20, help='min number of discordant read query bases flanking breakpoint on either side. (should not exceed read length) [20]')
    parser.add_argument('--split_weight', type=float, required=False, default=1, help='weight for split reads [1]')
    parser.add_argument('--disc_weight', type=float, required=False, default=1, help='weight for discordant paired-end reads [1]')
    parser.add_argument('-n', dest='num_samp', type=int, required=False, default=1000000, help='number of pairs to sample from BAM file for building insert size distribution [1000000]')
    parser.add_argument('-M', action='store_true', dest='legacy', required=False, help='split reads are flagged as secondary, not supplementary. For compatibility with legacy BWA-MEM "-M" flag')
    parser.add_argument('--debug', action='store_true', help='debugging verbosity')
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
