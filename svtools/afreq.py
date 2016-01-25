import argparse, sys
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant

class UpdateInfo(object):
    def __init__(self, vcf_stream):
        self.vcf_stream = vcf_stream

    def execute(self, output_handle=sys.stdout):
        in_header = True
        header = []
        breakend_dict = {} # cache to hold unmatched generic breakends for genotyping
        vcf = Vcf()
        vcf_out = output_handle

        # read input VCF
        for line in self.vcf_stream:
            if in_header:
                if line.startswith('##'):
                    header.append(line) 
                    continue
                elif line.startswith('#CHROM'):
                    v = line.rstrip().split('\t')
                    header.append('\t'.join(v[:8]))

                    in_header = False
                    vcf.add_header(header)
                    
                    vcf.add_info('AF', 'A', 'Float', 'Allele Frequency, for each ALT allele, in the same order as listed')
                    vcf.add_info('NSAMP', '1', 'Integer', 'Number of samples with non-reference genotypes')

                    # write header
                    vcf_out.write(vcf.get_header(include_samples=False))
                    vcf_out.write('\t' + '\t'.join(v[8:]) + '\n')
                continue

            v = line.rstrip().split('\t')
            var = Variant(v[:8], vcf)

            # extract genotypes from VCF
            num_alt = len(var.alt.split(','))
            alleles = [0] * (num_alt + 1)
            num_samp = 0

            for i in xrange(9,len(v)):
                gt_string = v[i].split(':')[0]

                if '.' in  gt_string:
                    continue
                gt = gt_string.split('/')
                if len(gt) == 1:
                    gt = gt_string.split('|')
                gt = map(int, gt)

                for i in xrange(len(gt)):
                    alleles[gt[i]] += 1

                # iterate the number of non-reference samples
                if sum(gt) > 0:
                    num_samp += 1

            allele_sum = float(sum(alleles))
            allele_freq = ['.'] * len(alleles)

            # populate AF
            if allele_sum > 0:
                for i in xrange(len(alleles)):
                    allele_freq[i] = alleles[i] / allele_sum
                var.info['AF'] = ','.join(map(str, ['%.4g' % a for a in allele_freq[1:]]))
            else:
                var.info['AF'] = ','.join(map(str, allele_freq[1:]))
            
            # populate NSAMP
            var.info['NSAMP'] = num_samp

            # after all samples have been processed, write
            vcf_out.write(var.get_var_string()
                          + '\t'
                          + '\t'.join(v[8:])
                          + '\n')
        vcf_out.close()

def description():
    return 'Add allele frequency information to a VCF file'

def add_arguments_to_parser(parser):
    parser.add_argument(metavar='vcf', dest='input_vcf', nargs='?', default=None, help='VCF input')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    handle = None
    if args.input_vcf == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            handle = sys.stdin
    else:
        handle = open(args.input_vcf, 'r')
    updater = UpdateInfo(handle)
    updater.execute()
    if handle != sys.stdin:
        handle.close()

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
