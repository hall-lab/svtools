import sys
import argparse
import svtyper.classic as typer

def description():
    return 'compute genotype of structural variants based on breakpoint depth'

def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--input_vcf', metavar='<VCF>', type=argparse.FileType('r'), default=None, help='VCF input (default: stdin)')
    parser.add_argument('-o', '--output_vcf', metavar='<VCF>', type=argparse.FileType('w'), default=sys.stdout, help='output VCF to write (default: stdout)')
    parser.add_argument('-B', '--bam', metavar='<BAM>', type=str, required=True, help='BAM or CRAM file')
    parser.add_argument('-T', '--ref_fasta', metavar='<FASTA>', type=str, required=False, default=None, help='Indexed reference FASTA file (recommended for reading CRAM files)')
    parser.add_argument('-l', '--lib_info', metavar='<JSON>', dest='lib_info_path', type=str, required=False, default=None, help='create/read JSON file of library information')
    parser.add_argument('-m', '--min_aligned', metavar='<INT>', type=int, required=False, default=20, help='minimum number of aligned bases to consider read as evidence [20]')
    parser.add_argument('-n', dest='num_samp', metavar='<INT>', type=int, required=False, default=1000000, help='number of pairs to sample from BAM file for building insert size distribution [1000000]')
    parser.add_argument('-q', '--sum_quals', action='store_true', required=False, help='add genotyping quality to existing QUAL (default: overwrite QUAL field)')
    parser.add_argument('--max_reads', metavar='<INT>', type=int, default=10000, required=False, help='maximum number of reads to assess at any variant (reduces processing time in high-depth regions, default: 10000)')
    parser.add_argument('--max_ci_dist', metavar='<INT>', type=int, default=1e10, required=False, help='maximum size of a confidence interval before 95%% CI is used intead (default: 1e10)')
    parser.add_argument('--split_weight', metavar='<FLOAT>', type=float, required=False, default=1, help='weight for split reads [1]')
    parser.add_argument('--disc_weight', metavar='<FLOAT>', type=float, required=False, default=1, help='weight for discordant paired-end reads [1]')
    parser.add_argument('-w', '--write_alignment', metavar='<BAM>', dest='alignment_outpath', type=str, required=False, default=None, help='write relevant reads to BAM file')
    parser.add_argument('--debug', action='store_true', help=argparse.SUPPRESS)
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    # if no input, check if part of pipe and if so, read stdin.
    if args.input_vcf == None:
        if not sys.stdin.isatty():
            args.input_vcf = sys.stdin

    typer.sv_genotype(args.bam,
                 args.input_vcf,
                 args.output_vcf,
                 args.min_aligned,
                 args.split_weight,
                 args.disc_weight,
                 args.num_samp,
                 args.lib_info_path,
                 args.debug,
                 args.alignment_outpath,
                 args.ref_fasta,
                 args.sum_quals,
                 args.max_reads,
                 args.max_ci_dist)

if __name__ == "__main__":
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
