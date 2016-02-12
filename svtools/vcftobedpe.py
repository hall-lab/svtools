import argparse
import sys
import time
import re

import svtools.vcf.file
import svtools.vcf.variant
import svtools.utils as su

def writeBND(prim, sec, v, bedpe_out):
    '''
    Parse out strand orientation from BND alt fields
    Simple mapping seems to be left-most strand corresponds to the direction of the square brackets. Right-most to side the reference base is on.

    For example:

    N[2:22222[ -> brackets pt left so - and N is on the left so plus +-
    ]2:22222]N -> brackets pt right so + and N is on the right so minus -+
    N]2:222222] -> brackets pt right so + and N is on the left so plus ++
    [2:222222[N -> brackets pt left so + and N is on the right so minus --
    '''
    primary = prim
    secondary = sec
    if prim is None:
        primary = sec 
    b1 = primary.pos
    o1 = '+'
    o2 = '+'
    sep = ']'
    if ']' not in primary.alt:
        sep = '['
        o2 = '-'
    if primary.alt.startswith('[') or primary.alt.startswith(']'):
            o1 = '-'
    r = re.compile(r'\%s(.+?)\%s' % (sep, sep))
    chrom2, b2 = r.findall(primary.alt)[0].split(':')
    b2 = int(b2)
    primary.set_info('END',b2)
    # XXX Colby mentioned that coordinates are 0 based sometimes
    # and not 0-based other times. This appears to be code pertaining to that
    # If reference is on the right of the brackets, o1 == '-' and coords adjusted
    if 'CIPOS' in primary.info:
        span = map(int, primary.info['CIPOS'].split(','))
        if o1 == '-':
            b1 -= 1
        s1 = b1 + span[0]
        e1 = b1 + span[1]
    else:
        if o1 == '-':
            b1 -= 1
        else:
            s1 = b1
            e1 = b1
    
    if 'CIEND' in primary.info:    
        span = map(int, primary.info['CIEND'].split(','))
        if o2 == '-':
            b2 -= 1
        s2 = b2 + span[0]
        e2 = b2 + span[1]
    else:
        if o2== '-':
            b2 -= 1
        else:
            s2 = b2
            e2 = b2

    ispan = s2 - e1
    ospan = e2 - s1
    chrom_A = primary.chrom
    chrom_B = chrom2
    
    # write bedpe
    #Swap fields for no primary present as we did calculation with secondary
    if prim is None:
        info_A = "MISSING"
        chrom_A = chrom2
        chrom_B = primary.chrom
        s1, s2 = s2, s1
        e1, e2 = e2, e1
        o1, o2 = o2, o1
    else:
        info_A = primary.get_info_string()    
    if sec is None:
        info_B = "MISSING"
    else:
        info_B = secondary.get_info_string()    
    bedpe_out.write('\t'.join(map(str,
                                  [chrom_A,
                                   max(s1,1) - 1,
                                   max(e1,0),
                                   chrom_B,
                                   max(s2,1) - 1,
                                   max(e2,0),
                                   primary.info['EVENT'],
                                   primary.original_qual,
                                   o1,
                                   o2,
                                   primary.info['SVTYPE'],
                                   primary.filter] + [info_A] + [info_B] + v[8:]
                                  )) + '\n')
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
            b1 = var.pos
            b2 = int(var.info['END'])
            name = v[2]
            if 'STRANDS' in var.info:
                strands = var.info['STRANDS']
                o1 = strands[0]
                o2 = strands[1]
            else:
                o1 = '+'
                o2 = '+'
            if 'CIPOS' in var.info:
                span = map(int, var.info['CIPOS'].split(','))
                s1 = b1 + span[0]
                e1 = b1 + span[1]
            else:
                e1 = b1
                s1 = b1
            if 'CIEND' in var.info:    
                span = map(int, var.info['CIEND'].split(','))
                s2 = b2 + span[0]
                e2 = b2 + span[1]
            else:
                e2 = b2
                s2 = b2    

            ispan = s2 - e1
            ospan = e2 - s1
            # write bedpe
            bedpe_out.write('\t'.join(map(str,
                                          [var.chrom,
                                           max(s1,0),
                                           max(e1,0),
                                           var.chrom,
                                           max(s2,0),
                                           max(e2,0),
                                           name,
                                           var.original_qual,
                                           o1,
                                           o2,
                                           var.info['SVTYPE'],
                                           var.filter] +
                                           [var.get_info_string()] + ['.'] + v[8:]
                                          )) + '\n')
        else:
            if 'SECONDARY' in var.info:
                if var.info['EVENT'] in bnds:
                    #primary
                    var1 = bnds[var.info['EVENT']]
                    writeBND(var1,var,v,bedpe_out)
                    del bnds[var.info['EVENT']]                              
                else:
                    sec_bnds.update({var.info['EVENT']:var})
            else: 
                bnds.update({var.info['EVENT']:var})
                continue
    intersected_keys = bnds.viewkeys() & sec_bnds.viewkeys()
    for key in intersected_keys:
       writeBND(bnds[key],sec_bnds[key],v,bedpe_out)
       del bnds[key] 
       del sec_bnds[key]
    if bnds is not None:
        for bnd in bnds:
            sys.stderr.write('Warning: missing secondary multiline variant at ID:' + bnds[bnd].info['EVENT'] + '\n')
            writeBND(bnds[bnd],None,v,bedpe_out)
    if sec_bnds is not None:
        for bnd in sec_bnds:
            sys.stderr.write('Warning: missing primary multiline variant at ID:' + sec_bnds[bnd].info['EVENT'] + '\n')
            writeBND(None,sec_bnds[bnd],v,bedpe_out)
            
    # close the files
    bedpe_out.close()
    vcf_file.close()

    return

def description():
    return 'Convert a VCF file to a BEDPE file'

def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--input', default=None, help='VCF input (default: stdin)')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help='Output BEDPE to write (default: stdout)')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
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
