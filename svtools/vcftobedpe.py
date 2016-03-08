import argparse
import sys
import time
import re

import svtools.vcf.file
import svtools.vcf.variant
import svtools.utils as su

def convert_simple(vcf_line_array, variant):
    '''
    Convert simple SVs to BEDPE
    Precise breakpoints are specified by 0-length intervals.
    For simple SVs, breakpoint 1 is the VCF coordinate and breakpoint 2 is the END coordinate
    '''
    start1 = end1 = variant.pos
    start2 = end2 = int(variant.info['END'])
    orientation1 = orientation2 = '+'
    if 'STRANDS' in variant.info:
        strands = variant.info['STRANDS']
        orientation1, orientation2 = strands[:2]
    if 'CIPOS' in variant.info:
        span = map(int, variant.info['CIPOS'].split(','))
        start1 += span[0]
        end1 += span[1]
    if 'CIEND' in variant.info:
        span = map(int, variant.info['CIEND'].split(','))
        start2 += span[0]
        start2 += span[1]
    return '\t'.join(map(str,
        [variant.chrom,
        max(start1, 0),
        max(end1, 0),
        variant.chrom,
        max(start2, 0),
        max(end2, 0),
        variant.var_id,
        variant.original_qual,
        orientation1,
        orientation2,
        variant.info['SVTYPE'],
        variant.filter] +
        [variant.get_info_string()] + ['.'] + vcf_line_array[8:]
        ))

def convert_breakend(vcf_line_array, primary_variant, secondary_variant):
    '''
    Parse out strand orientation from BND alt fields
    Simple mapping seems to be left-most strand corresponds to the direction of the square brackets. Right-most to side the reference base is on.

    For example:

    N[2:22222[ -> brackets pt left so - and N is on the left so plus +-
    ]2:22222]N -> brackets pt right so + and N is on the right so minus -+
    N]2:222222] -> brackets pt right so + and N is on the left so plus ++
    [2:222222[N -> brackets pt left so + and N is on the right so minus --
    '''
    variant = primary_variant
    if primary_variant is None:
        variant = secondary_variant
    chrom1 = variant.chrom
    breakpoint1 = variant.pos
    orientation1 = '+'
    orientation2 = '+'
    # NOTE The below is ugly but intended to match things like [2:222[ and capture the brackets
    r = re.compile(r'([][])(.+?)([][])')
    sep1, region, sep2 = r.findall(variant.alt)[0]
    
    assert sep1 == sep2
    assert sep1 != None
    
    chrom2, breakpoint2 = region.split(':')
    breakpoint2 = int(breakpoint2)

    if variant.alt.startswith(sep1):
        orientation1 = '-'
        breakpoint1 -= 1
    if sep1 == '[':
        orientation2 = '-'
        breakpoint2 -= 1

    start1 = end1 = breakpoint1
    start2 = end2 = breakpoint2
    variant.set_info('END', breakpoint2)

    if 'CIPOS' in variant.info:
        span = map(int, variant.info['CIPOS'].split(','))
        start1 += span[0]
        end1 += span[1]
    
    if 'CIEND' in variant.info:    
        span = map(int, variant.info['CIEND'].split(','))
        start2 += span[0]
        end2 += span[1]

    chrom_A = variant.chrom
    chrom_B = chrom2
    
    # write bedpe
    #Swap fields for no primary present as we did calculation with secondary
    if primary_variant is None:
        info_A = "MISSING"
        chrom_A, chrom_B = chrom_B, chrom_A
        start1, start2 = start2, start1
        end1, end2 = end2, end1
        orientation1, orientation2 = orientation2, orientation1
    else:
        info_A = primary_variant.get_info_string()    
    if secondary_variant is None:
        info_B = "MISSING"
    else:
        info_B = secondary_variant.get_info_string()    
    return '\t'.join(map(str,
                                  [chrom_A,
                                   max(start1,0),
                                   max(end1,0),
                                   chrom_B,
                                   max(start2,0),
                                   max(end2,0),
                                   variant.info['EVENT'],
                                   variant.original_qual,
                                   orientation1,
                                   orientation2,
                                   variant.info['SVTYPE'],
                                   variant.filter] + [info_A] + [info_B] + vcf_line_array[8:]
                                  ))

# primary function
def vcfToBedpe(vcf_file, bedpe_out):
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
                continue
            else:
                # print header
                in_header = False
                vcf.add_header(header)
                if "SVTYPE" in [info.id for info in vcf.info_list]:
                   vcf.add_info_after("SVTYPE", "POS", 1, 'Integer', 'Position of the variant described in this record')
                header=vcf.get_header()
                bedpe_out.write(header[:header.rfind('\n')] + '\n')                
                if len(sample_list) > 0:
                    bedpe_out.write('\t'.join(['#CHROM_A',
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
                                               'INFO_A','INFO_B',
                                               'FORMAT','\t'.join(map(str,sample_list))] 
                                             ) + '\n')
                else:
                    bedpe_out.write('\t'.join(['#CHROM_A',
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
                                               'INFO_A','INFO_B']
                                              ) + '\n')

        v = line.rstrip().split('\t')
        var = svtools.vcf.variant.Variant(v, vcf)
        var.set_info("POS", var.pos)
        if var.info['SVTYPE'] != 'BND':
            bedpe_out.write(convert_simple(v, var))
            bedpe_out.write('\n')
        else:
            if 'SECONDARY' in var.info:
                if var.info['EVENT'] in bnds:
                    #primary
                    var1 = bnds[var.info['EVENT']]
                    bedpe_out.write(convert_breakend(v, var1, var))
                    bedpe_out.write('\n')
                    del bnds[var.info['EVENT']]                              
                else:
                    sec_bnds.update({var.info['EVENT']:var})
            else: 
                bnds.update({var.info['EVENT']:var})
                continue
    intersected_keys = bnds.viewkeys() & sec_bnds.viewkeys()
    for key in intersected_keys:
        bedpe_out.write(convert_breakend(v, bnds[key], sec_bnds[key]))
        bedpe_out.write('\n')
        del bnds[key] 
        del sec_bnds[key]
    if bnds is not None:
        for bnd in bnds:
            sys.stderr.write('Warning: missing secondary multiline variant at ID:' + bnds[bnd].info['EVENT'] + '\n')
            bedpe_out.write(convert_breakend(v, bnds[bnd], None))
            bedpe_out.write('\n')
    if sec_bnds is not None:
        for bnd in sec_bnds:
            sys.stderr.write('Warning: missing primary multiline variant at ID:' + sec_bnds[bnd].info['EVENT'] + '\n')
            bedpe_out.write(convert_breakend(v, None, sec_bnds[bnd]))
            bedpe_out.write('\n')
            
    # close the files
    bedpe_out.close()
    vcf_file.close()

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
