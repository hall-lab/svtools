import argparse, sys, re
import math, time
import copy
from argparse import RawTextHelpFormatter
from svtools.bedpe import Bedpe
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
import svtools.utils as su

def adjust_by_cipos(bedpe):
    '''
    Undo adjustments to BEDPE coordinates for the VCF POS based on confidence intervals
    '''
    position = bedpe.s1
    if bedpe.o1 == '-' and bedpe.svtype == 'BND':
        position += 1   #undo left adjust based on strandedness
    if 'CIPOS=' in bedpe.info1:
        cipos = re.split('=|,', ''.join(filter(lambda x: 'CIPOS=' in x, bedpe.info1.split(';'))))
        position -= int(cipos[1])
    return position

def adjust_by_ciend(bedpe):
    '''
    Undo adjustments to bEDPE coordinates for the VCF END field based on confidence intervals
    '''
    end = bedpe.s2
    if bedpe.o2 == '-' and bedpe.svtype == 'BND':
        end += 1

    if 'CIEND=' in bedpe.info1:     
        ciend = re.split('=|,', ''.join(filter(lambda x: 'CIEND=' in x, bedpe.info1.split(';'))))
        end -= int(ciend[1])
    return end

# primary function
def bedpeToVcf(bedpe_file, vcf_out):
    myvcf = Vcf()
    in_header = True
    # parse the bedpe data
    header = list()
    for line in bedpe_file:
        if in_header:
            if line[0:2] == '##':
                header.append(line)
                continue
            elif line[0] == '#' and line[1] != '#':    
                sample_list_str = line.rstrip().split('\t', 14)[-1]
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
        b1 = adjust_by_cipos(bedpe)
        b2 = adjust_by_ciend(bedpe)
        if bedpe.svtype == 'BND':
            bedpe1_list = [
                    bedpe.c1, 
                    b1,
                    bedpe.name + '_1', #ID
                    'N',
                    '<' + str(bedpe.svtype) + '>', #ALT
                    bedpe.score,
                    bedpe.filter,
                    bedpe.info1
                    ]
            bedpe1_list.extend(bedpe.misc)
            var1 = Variant(bedpe1_list, myvcf)
            if bedpe.o1 == '+':
                if bedpe.o2 == '-':
                    var1.alt = '%s[%s:%s[' % (var1.ref, bedpe.c2, b2)
                elif bedpe.o2 == '+':
                    var1.alt = '%s]%s:%s]' % (var1.ref, bedpe.c2, b2)
            elif bedpe.o1 == '-':
                if bedpe.o2 == '+':
                    var1.alt = ']%s:%s]%s' % (bedpe.c2, b2, var1.ref)
                elif bedpe.o2 == '-':
                    var1.alt = '[%s:%s[%s' % (bedpe.c2, b2, var1.ref)
            bedpe_info = bedpe.info1
            info = copy.deepcopy(bedpe.info1)
            strands = re.split('=|:',''.join(filter(lambda x: 'STRANDS=' in x, bedpe_info.split(";"))))
            strands_str = str(strands[0]) + '=' + str(strands[1][::-1]) + ':' + str(strands[2])
            info = info.replace(''.join(filter(lambda x: 'STRANDS=' in x, bedpe.info1.split(";"))), strands_str)
            #add the cipos ciend,cipos95 and ciend95 variables
            info = info.replace(''.join(filter(lambda x: 'CIPOS=' in x, bedpe_info.split(";"))),'CIPOS='+ re.split('=',''.join(filter(lambda x: 'CIEND=' in x, bedpe_info.split(";"))))[1])            
            info = info.replace(''.join(filter(lambda x: 'CIEND='  in x, bedpe_info.split(";"))),'CIEND='+ re.split('=',''.join(filter(lambda x: 'CIPOS=' in x, bedpe_info.split(";"))))[1])
            info = info.replace(''.join(filter(lambda x: 'CIPOS95=' in x, bedpe_info.split(";"))),'CIPOS95='+ re.split('=',''.join(filter(lambda x: 'CIEND95=' in x, bedpe_info.split(";"))))[1])
            info = info.replace(''.join(filter(lambda x: 'CIEND95=' in x, bedpe_info.split(";"))),'CIEND95='+ re.split('=',''.join(filter(lambda x: 'CIPOS95=' in x, bedpe_info.split(";"))))[1])
            #Change MATEID
            info = info.replace(''.join(filter(lambda x: 'MATEID=' in x, bedpe_info.split(";"))),'MATEID=' + bedpe.name + '_2')
            #ADD IDENTIFIER FOR SECONDARY BREAKEND MATE
            info = info.replace(''.join(filter(lambda x: 'EVENT=' in x, bedpe_info.split(";"))),''.join(filter(lambda x: 'EVENT=' in x, bedpe_info.split(";"))) + ';SECONDARY;')

            bedpe2_list = [
                    bedpe.c2,  #chrom1
                    b2,
                    bedpe.name + '_2', #ID
                    'N',
                    '<' + str(bedpe.svtype) + '>', #ALT
                    bedpe.score,
                    bedpe.filter,
                    info
                    ]
            bedpe2_list.extend(bedpe.misc)

            var2 = Variant(bedpe2_list, myvcf)
            # add the strands field. For variant 2 must switch the order
            if bedpe.o2 == '+':
                if bedpe.o1 == '-':
                    var2.alt = '%s[%s:%s[' % (var2.ref, bedpe.c1, b1)
                elif bedpe.o1 == '+':
                    var2.alt = '%s]%s:%s]' % (var2.ref, bedpe.c1, b1)
            elif bedpe.o2 == '-':
                if bedpe.o1 == '+':
                    var2.alt = ']%s:%s]%s' % (bedpe.c1, b1, var2.ref)
                elif bedpe.o1 == '-':
                    var2.alt = '[%s:%s[%s' % (bedpe.c1, b1, var2.ref)
            if bedpe.malformedFlag == 0:
                vcf_out.write(var1.get_var_string() + '\n')
                vcf_out.write(var2.get_var_string() + '\n')
            elif bedpe.malformedFlag == 1:
                vcf_out.write(var2.get_var_string() + '\n')
            elif bedpe.malformedFlag == 2:
                vcf_out.write(var1.get_var_string() + '\n')
        else:
            # set VCF info elements for simple events
            bedpe_list = [
                    bedpe.c1,  #chrom1
                    b1,
                    bedpe.name, #ID
                    'N',
                    '<' + str(bedpe.svtype) + '>', #ALT
                    bedpe.score,
                    bedpe.filter,
                    bedpe.info1
                    ]
            bedpe_list.extend(bedpe.misc)

            var = Variant(bedpe_list, myvcf)
            # write the record to the VCF output file
            vcf_out.write(var.get_var_string() + '\n')

    # close the VCF output file
    vcf_out.close()
    
    return

def description():
    return 'convert a BEDPE file to VCF'

def epilog():
    return 'The input BEDPE file can be gzipped if it is specified explicitly.'

def add_arguments_to_parser(parser):
    parser.add_argument('-b', '--bedpe', metavar='<BEDPE>', default=None, help='BEDPE input (default: stdin)')
    parser.add_argument('-o', '--output', metavar='<VCF>', type=argparse.FileType('w'), default=sys.stdout, help='Output VCF to write (default: stdout)')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.bedpe) as stream:
        bedpeToVcf(stream, args.output)

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
