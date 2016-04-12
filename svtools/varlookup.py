import sys
import gzip
from operator import itemgetter
import argparse
from svtools.vcf.file import Vcf
from svtools.bedpe import Bedpe

def get_var_string(bedpe, cohort_name):
    if len(bedpe.cohort_vars) > 0:
        bedpe.set_info(cohort_name + '_AF', ','.join([value for (key, value) in sorted(bedpe.cohort_vars.items(), key=itemgetter(1),reverse=True)]))
        bedpe.set_info(cohort_name + '_VarID', ','.join([key for (key, value) in sorted(bedpe.cohort_vars.items(), key=itemgetter(1),reverse=True)]))
    else:
        bedpe.set_info(cohort_name + '_AF', str(0))
        bedpe.set_info(cohort_name + '_VarID', 'NONE')
    return str(bedpe)

def add(a_bedpe, b_bedpe, max_distance):
    if a_bedpe.svtype ==  b_bedpe.svtype:
        if (a_bedpe.o1 != b_bedpe.o1
            or a_bedpe.o2 != b_bedpe.o2):
            return False

        if (a_bedpe.c1 != b_bedpe.c1
            or a_bedpe.s1 - max_distance > b_bedpe.e1
            or a_bedpe.e1 + max_distance < b_bedpe.s1):
            return False

        if (a_bedpe.c2 != b_bedpe.c2
            or a_bedpe.s2 - max_distance > b_bedpe.e2
            or a_bedpe.e2 + max_distance < b_bedpe.s2):
            return False
        else:
            a_bedpe.cohort_vars[b_bedpe.name] = b_bedpe.af
            return True
    else:
        return False

def varLookup(aFile, bFile, bedpe_out, max_distance, pass_prefix, cohort_name):
    # FIXME The following code is heavily duplicated with vcftobedpe and bedpetovcf. Harmonize!!!
    bList = list()
    headerObj=Vcf() #co-opt the VCF header object
    if cohort_name is None:
        cohort_name=str(str(bFile).split('/')[-1])
        
    if bFile == "stdin":
        bData = sys.stdin
    elif bFile.endswith('.gz'):
        bData = gzip.open(bFile, 'rb')
    else:
        bData = open(bFile, 'r')
    for bLine in bData:
        if bLine.startswith(pass_prefix):
            continue
        bentry = Bedpe(bLine.rstrip().split('\t'))
        if bentry.af is None:
            sys.stderr.write('No allele frequency for variant found in -b file. This tool requires allele frequency information to function. Please add with svtools afreq and rerun\n')
            sys.exit(1)
        bList.append(bentry)
    
    if aFile == "stdin":
        aData = sys.stdin
    elif aFile.endswith('.gz'):
        aData = gzip.open(aFile, 'rb')
    else:
        aData = open(aFile, 'r')
    in_header=True    
    header_lines = []
    sample_list = None
    for aLine in aData:
        if pass_prefix is not None and aLine.startswith(pass_prefix):
            if aLine[0] == '#' and aLine[1] != '#':
                sample_list = aLine.rstrip().split('\t', 14)[-1]
            else:
                header_lines.append(aLine)
            continue
        else:
            if in_header == True:
                headerObj.add_header(header_lines)
                headerObj.add_info(cohort_name + '_AF', '.', 'Float', 'Allele frequency(ies) for matching variants found in the ' + cohort_name + ' vcf' + ' (' + str(str(bFile).split('/')[-1]) + ')' )
                headerObj.add_info(cohort_name + '_VarID', '.', 'Integer', 'List of Variant ID(s) for matching variants found in the ' + cohort_name + ' vcf' + ' (' + str(str(bFile).split('/')[-1]) + ')' )

                header = headerObj.get_header()
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
                                               sample_list]
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
                in_header=False
            a = Bedpe(aLine.rstrip().split('\t'))
            if a.af is None:
                sys.stderr.write('No allele frequency for variant found in -a file. This tool requires allele frequency information to function. Please add with svtools afreq and rerun\n')
                sys.exit(1)
            for b in bList:
                add(a,b,max_distance)
            bedpe_out.write(get_var_string(a, cohort_name) + '\n')

def description():
    return 'look for variants common between two BEDPE files'

def add_arguments_to_parser(parser):
    parser.add_argument('-d', '--distance', metavar='<INT>', type=int, dest='max_distance', default=50, help='max separation distance (bp) of adjacent loci between bedpe files [50]')
    parser.add_argument("-a", "--aFile", dest="aFile", metavar='<BEDPE>', help="pruned, merged BEDPE (A file) or standard input (-a stdin).")
    parser.add_argument("-b", "--bFile", dest="bFile", metavar='<BEDPE>', help="pruned merged BEDPE (B file) (-b stdin). For pruning use svtools prune")
    parser.add_argument("-c", "--cohort", dest='cohort_name', metavar='<STRING>', default=None, help="cohort name to add information of matching variants (default:bFile)")                    
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<BEDPE>', default=sys.stdout, help='output BEDPE to write (default: stdout)')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    pass_prefix = "#"
    if args.aFile == None:
        if sys.stdin.isatty():
            sys.stderr.write('Please stream in input to this command or specify the file to read\n')
            sys.exit(1)
        else:
            args.aFile = sys.stdin

    try:
        varLookup(args.aFile, args.bFile, args.output, args.max_distance, pass_prefix, args.cohort_name)
    except IOError as err:
        sys.stderr.write("IOError " + str(err) + "\n");

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
