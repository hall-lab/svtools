import argparse, sys
from svtools.bedpe import Bedpe
from svtools.vcf.file import Vcf
from svtools.bedpetovcfconverter import BedpeToVcfConverter
import svtools.utils as su


# primary function
def bedpeToVcf(bedpe_file, vcf_out):
    myvcf = Vcf()
    converter = BedpeToVcfConverter(myvcf)
    in_header = True
    # parse the bedpe data
    header = list()
    for line in bedpe_file:
        if in_header:
            if line[0:2] == '##':
                header.append(line)
                continue
            elif line[0] == '#' and line[1] != '#':
                sample_list_str = line.rstrip().split('\t', 20)[-1]
                header.append('\t'.join([
                                    '#CHROM',
                                    'POS',
                                    'ID',
                                    'REF',
                                    'ALT',
                                    'QUAL',
                                    'FILTER',
                                    'INFO',
                                    sample_list_str
                                    ] ))
                continue
            else:
                in_header = False
                myvcf.add_header(header)
                myvcf.file_format='VCFv4.2'
                vcf_out.write(myvcf.get_header() + '\n')
        #
        bedpe = Bedpe(line.rstrip().split('\t'))
        variants = converter.convert(bedpe)
        for v in variants:
            vcf_out.write(v.get_var_string() + '\n')

    # close the VCF output file and header if no variants found
    if in_header == True:
        myvcf.add_header(header)
        myvcf.file_format='VCFv4.2'
        vcf_out.write(myvcf.get_header() + '\n')
    vcf_out.close()

    return

def description():
    return 'convert a BEDPE file to VCF'

def epilog():
    return 'The input BEDPE file can be gzipped if it is specified explicitly.'

def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--input', metavar='<BEDPE>', default=None, help='BEDPE input (default: stdin)')
    parser.add_argument('-o', '--output', metavar='<VCF>', type=argparse.FileType('w'), default=sys.stdout, help='Output VCF to write (default: stdout)')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.input) as stream:
        bedpeToVcf(stream, args.output)

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
