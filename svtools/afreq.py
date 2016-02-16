import argparse, sys
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
import svtools.utils as su

class UpdateInfo(object):
    def __init__(self, vcf_stream):
        self.vcf_stream = vcf_stream

    def calc_msq(self, var):
            # Below is what was in vcfpaste, but what if multiple ALTs?
            # Do we only expect 1 ALT per line?
            # NOTE SQ is defined as: 'Phred-scaled probability that this site is variant (non-reference in this sample'
            # Likely want average sample quality across all non-0/0 genotypes rather than just those containing 1
            gt = [var.genotype(s).get_format('GT') for s in var.sample_list]
            positive_gt = filter(lambda x: x == '0/1' or x == '1/1', gt)
            num_pos = len(positive_gt)
            sum_sq = 0.0
            try:
                sum_sq += sum([float(var.genotype(s).get_format('SQ')) for s in var.sample_list \
                if var.genotype(s).get_format('GT') == '1/1' or var.genotype(s).get_format('GT') == '0/1'])
            except ValueError:
                sum_sq += 0
            if num_pos > 0:
                msq = '%0.2f' % (sum_sq / num_pos)
            else:
                msq = '.'
            return msq

    def execute(self, output_handle=sys.stdout):
        in_header = True
        header = []
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
                    header.append('\t'.join(v))

                    in_header = False
                    vcf.add_header(header)
                    
                    vcf.add_info('AF', 'A', 'Float', 'Allele Frequency, for each ALT allele, in the same order as listed')
                    vcf.add_info('NSAMP', '1', 'Integer', 'Number of samples with non-reference genotypes')
                    vcf.add_info('MSQ', '1', 'Float', 'Mean sample quality of positively genotyped samples')

                    # write header
                    vcf_out.write(vcf.get_header() + '\n')
                    #vcf_out.write('\t' + '\t'.join(v[8:]) + '\n')
                continue

            v = line.rstrip().split('\t')
            var = Variant(v, vcf, fixed_genotypes=True)

            # extract genotypes from VCF
            num_alt = len(var.alt.split(','))
            alleles = [0] * (num_alt + 1)
            num_samp = 0

            gt = [var.genotype(s).get_format('GT') for s in var.sample_list]
            for gt_string in gt:

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
            var.info['MSQ'] = self.calc_msq(var)

            # after all samples have been processed, write
            vcf_out.write(var.get_var_string(use_cached_gt_string=True) + '\n')
        vcf_out.close()

def description():
    return 'add allele frequency information to a VCF file'

def epilog():
    return 'Specify the path to an (optionally) bgzipped VCF. If no file is specified then input is read from stdin.'


def add_arguments_to_parser(parser):
    parser.add_argument(metavar='<VCF>', dest='input_vcf', nargs='?', default=None, help='VCF input')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description(), epilog=epilog())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.input_vcf) as input_stream:
        updater = UpdateInfo(input_stream)
        updater.execute()

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
