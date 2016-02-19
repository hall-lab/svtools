import argparse, sys, re
import math, time
import copy
from argparse import RawTextHelpFormatter
from svtools.bedpe import Bedpe
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
import svtools.utils as su

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
        if bedpe.svtype == 'BND':
            bedpe1_list = [
                    bedpe.c1, 
                    bedpe.b1,
                    bedpe.name + '_1', #ID
                    'N',
                    '<' + str(bedpe.svtype) + '>', #ALT
                    bedpe.score,
                    bedpe.filter
                    ]
            bedpe1_list.extend(bedpe.misc)
            var1 = Variant(bedpe1_list, myvcf)
            if bedpe.o1 == '+':
                if bedpe.o2 == '-':
                    var1.alt = '%s[%s:%s[' % (var1.ref, bedpe.c2, bedpe.b2)
                elif bedpe.o2 == '+':
                    var1.alt = '%s]%s:%s]' % (var1.ref, bedpe.c2, bedpe.b2)
            elif bedpe.o1 == '-':
                if bedpe.o2 == '+':
                    var1.alt = ']%s:%s]%s' % (bedpe.c2, bedpe.b2, var1.ref)
                elif bedpe.o2 == '-':
                    var1.alt = '[%s:%s[%s' % (bedpe.c2, bedpe.b2, var1.ref)
            misc = copy.deepcopy(bedpe.misc)
            strands = re.split('=|:',''.join(filter(lambda x: 'STRANDS=' in x, bedpe.misc[0].split(";"))))
            strands_str = str(strands[0]) + '=' + str(strands[1][::-1]) + ':' + str(strands[2])
            misc[0]=misc[0].replace(''.join(filter(lambda x: 'STRANDS=' in x, bedpe.misc[0].split(";"))), strands_str)
            #add the cipos ciend,cipos95 and ciend95 variables
            misc[0]=misc[0].replace(''.join(filter(lambda x: 'CIPOS=' in x, bedpe.misc[0].split(";"))),'CIPOS='+ re.split('=',''.join(filter(lambda x: 'CIEND=' in x, bedpe.misc[0].split(";"))))[1])            
            misc[0]=misc[0].replace(''.join(filter(lambda x: 'CIEND='  in x, bedpe.misc[0].split(";"))),'CIEND='+ re.split('=',''.join(filter(lambda x: 'CIPOS=' in x, bedpe.misc[0].split(";"))))[1])
            misc[0]=misc[0].replace(''.join(filter(lambda x: 'CIPOS95=' in x, bedpe.misc[0].split(";"))),'CIPOS95='+ re.split('=',''.join(filter(lambda x: 'CIEND95=' in x, bedpe.misc[0].split(";"))))[1])
            misc[0]=misc[0].replace(''.join(filter(lambda x: 'CIEND95=' in x, bedpe.misc[0].split(";"))),'CIEND95='+ re.split('=',''.join(filter(lambda x: 'CIPOS95=' in x, bedpe.misc[0].split(";"))))[1])
            #Change MATEID
            misc[0]= misc[0].replace(''.join(filter(lambda x: 'MATEID=' in x, bedpe.misc[0].split(";"))),'MATEID=' + bedpe.name + '_2')
            #ADD IDENTIFIER FOR SECONDARY BREAKEND MATE
            misc[0]=misc[0].replace(''.join(filter(lambda x: 'EVENT=' in x, bedpe.misc[0].split(";"))),''.join(filter(lambda x: 'EVENT=' in x, bedpe.misc[0].split(";"))) + ';SECONDARY;')

            bedpe2_list = [
                    bedpe.c2,  #chrom1
                    bedpe.b2,
                    bedpe.name + '_2', #ID
                    'N',
                    '<' + str(bedpe.svtype) + '>', #ALT
                    bedpe.score,
                    bedpe.filter
                    ]
            bedpe2_list.extend(misc)

            var2 = Variant(bedpe2_list, myvcf)
            # add the strands field. For variant 2 must switch the order
            if bedpe.o2 == '+':
                if bedpe.o1 == '-':
                    var2.alt = '%s[%s:%s[' % (var2.ref, bedpe.c1, bedpe.b1)
                elif bedpe.o1 == '+':
                    var2.alt = '%s]%s:%s]' % (var2.ref, bedpe.c1, bedpe.b1)
            elif bedpe.o2 == '-':
                if bedpe.o1 == '+':
                    var2.alt = ']%s:%s]%s' % (bedpe.c1, bedpe.b1, var2.ref)
                elif bedpe.o1 == '-':
                    var2.alt = '[%s:%s[%s' % (bedpe.c1, bedpe.b1, var2.ref)
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
                    bedpe.b1,
                    bedpe.name, #ID
                    'N',
                    '<' + str(bedpe.svtype) + '>', #ALT
                    bedpe.score,
                    bedpe.filter
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
