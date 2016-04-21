import argparse, sys
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
import svtools.utils as su

class UpdateInfo(object):
    def __init__(self, vcf_stream):
        self.vcf_stream = vcf_stream

    @staticmethod
    def numeric_alleles(gt_string):
        '''
        Convert a VCF allele string into integers
        '''
        gt = gt_string.split('/')
        if len(gt) == 1:
            gt = gt_string.split('|')
        return map(int, gt)

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
            var = Variant(v, vcf)

            # extract genotypes from VCF
            num_alt = len(var.alt.split(','))
            alleles = [0] * (num_alt + 1)
            num_samp = 0
            sum_sq = 0.0

            for gt in var.genotypes():
                gt_string = gt.get_format('GT')

                if '.' not in gt_string:
                    indexes = self.numeric_alleles(gt_string)

                    for i in indexes:
                        alleles[i] += 1

                    # iterate the number of non-reference samples
                    if sum(indexes) > 0:
                        num_samp += 1
                        try:
                            sum_sq += float(gt.get_format('SQ'))
                        except KeyError:
                            pass

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
            if num_samp > 0:
                msq = '%0.2f' % (sum_sq / num_samp)
            else:
                msq = '.'
            var.info['MSQ'] = msq

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
