import argparse
import sys
import time

import svtools.vcf.file
import svtools.vcf.variant
import svtools.utils as su
from svtools.vcftobedpeconverter import VcfToBedpeConverter

def vcfToBedpe(vcf_file, bedpe_out):
    converter = VcfToBedpeConverter()
    vcf = svtools.vcf.file.Vcf()
    in_header = True
    header = []
    sample_list = []
    bnds = dict()
    sec_bnds = dict()
    v = []
    for line in vcf_file:
        if in_header:
            if line[0:2] == '##':
                if line.split('=')[0] == '##fileformat':
                    line = '##fileformat=' + "BEDPE" + '\n'
                if line.split('=')[0] == '##fileDate':
                    line = '##fileDate=' + time.strftime('%Y%m%d') + '\n'
                header.append(line)
                continue
            elif line[0] == '#' and line[1] != '#':
                sample_list = line.rstrip().split('\t')[9:]
                header.append(line)
                continue
            else:
                # print header
                in_header = False
                vcf.add_header(header)
                if "SVTYPE" in [info.id for info in vcf.info_list]:
                   vcf.add_info_after("SVTYPE", "POS", 1, 'Integer', 'Position of the variant described in this record')
                header=vcf.get_header()
                bedpe_out.write(header[:header.rfind('\n')] + '\n')
                final_header_line = ['#CHROM_A',
                        'START_A',
                        'END_A',
                        'CHROM_B',
                        'START_B',
                        'END_B',
                        'ID',
                        'QUAL',
                        'STRAND_A',
                        'STRAND_B',
                        'TYPE',
                        'FILTER',
                        'NAME_A',
                        'REF_A',
                        'ALT_A',
                        'NAME_B',
                        'REF_B',
                        'ALT_B',
                        'INFO_A',
                        'INFO_B']

                if len(sample_list) > 0:
                    bedpe_out.write('\t'.join(final_header_line + ['FORMAT','\t'.join(map(str,sample_list))]) + '\n')
                else:
                    bedpe_out.write('\t'.join(final_header_line) + '\n')

        v = line.rstrip().split('\t')
        var = svtools.vcf.variant.Variant(v, vcf)
        var.set_info("POS", var.pos)
        # If there is no MATEID then assume this is a single-ended BND and simply output
        if var.info['SVTYPE'] != 'BND' or 'MATEID' not in var.info:
            bedpe_out.write(str(converter.convert(var)) + '\n')
        else:
            mate_id = var.info['MATEID']
            if 'SECONDARY' in var.info:
                if mate_id in bnds:
                    #primary
                    var1 = bnds[mate_id]
                    bedpe_out.write(str(converter.convert(var1, var)) + '\n')
                    del bnds[mate_id]
                else:
                    sec_bnds.update({var.var_id:var})
            else:
                if mate_id in sec_bnds:
                    var2 = sec_bnds[mate_id]
                    bedpe_out.write(str(converter.convert(var, var2)) + '\n')
                    del sec_bnds[mate_id]
                else:
                    bnds.update({var.var_id:var})
    if bnds is not None:
        for bnd in bnds:
            sys.stderr.write('Warning: missing secondary multiline variant at ID:' + bnd + '\n')
            bedpe_out.write(str(converter.convert(bnds[bnd], None)) + '\n')
    if sec_bnds is not None:
        for bnd in sec_bnds:
            sys.stderr.write('Warning: missing primary multiline variant at ID:' + bnd + '\n')
            bedpe_out.write(str(converter.convert(None, sec_bnds[bnd])) + '\n')

    # close the files
    bedpe_out.close()
    return

def description():
    return 'convert a VCF file to a BEDPE file'

def epilog():
    return 'The input VCF file can be gzipped if it is specified explicitly.'

def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--input', metavar='<VCF>', default=None, help='VCF input (default: stdin)')
    parser.add_argument('-o', '--output', metavar='<BEDPE>', type=argparse.FileType('w'), default=sys.stdout, help='output BEDPE to write (default: stdout)')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description(), epilog=epilog())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.input) as stream:
        return vcfToBedpe(stream, args.output)

# initialize the script
if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
