import sys
import string
import gzip
import time
from operator import itemgetter
from optparse import OptionParser
import argparse, sys , re
from svtools.vcf.file import Vcf

class Bedpe(object):
    def __init__(self, line):
        v = line.strip().split("\t")
        self.chrom_a = v[0]
        self.start_a = int(v[1])
        self.end_a = int(v[2])
        self.chrom_b = v[3]
        self.start_b = int(v[4])
        self.end_b = int(v[5])
        self.id = v[6]
        self.qual = v[7]
        self.info = v[12]
        self.af=filter(lambda x: x if x.startswith('AF=') else None, self.info.split(";"))
        #print self.af
        self.cohort_vars = dict()
        if self.af is not None and len(self.af) == 1:
            self.af=''.join(self.af).replace('AF=','')
        else: 
            print "No allele frequency for variants found. Input a valid file"
            sys.exit(0)
        
        self.sv_event=filter(lambda x: x if x.startswith('SVTYPE=') else None, self.info.split(";"))
        if self.sv_event is not None and len(self.sv_event) == 1:
            self.sv_event=''.join(self.sv_event).replace('SVTYPE=','')
        else:
            print "No svtype field found. Input a valid file"
            sys.exit(0)
        self.filter = v[11]
        try:
            self.strand_a = v[8]
            self.strand_b = v[9]
        except IndexError:
            self.strand_a = ''
            self.strand_b = ''
        self.vec = '\t'.join(v[13:])
    def get_var_string(self,cohort_name):
        if len(self.cohort_vars)>0:
            self.info = self.info + ';' + cohort_name + '_AF=' + ','.join([value for (key, value) in sorted(self.cohort_vars.items(), key=itemgetter(1),reverse=True)])
            self.info = self.info + ';' + cohort_name + '_VarID=' + ','.join([key for (key, value) in sorted(self.cohort_vars.items(), key=itemgetter(1),reverse=True)])
        else:
            self.info = self.info + ';' + cohort_name + '_AF=' + str(0)
            self.info = self.info + ';' + cohort_name + '_VarID=' + 'NONE'
        return '\t'.join([self.chrom_a,str(self.start_a),str(self.end_a),\
                        self.chrom_b,str(self.start_b),str(self.end_b),\
                        self.id,self.qual,self.strand_a,self.strand_b,\
                        self.sv_event,self.filter,self.info,self.vec]) + '\n'

def add(a_bedpe,b_bedpe,max_distance):
    if a_bedpe.sv_event ==  b_bedpe.sv_event:
        if (a_bedpe.strand_a != b_bedpe.strand_a
            or a_bedpe.strand_b != b_bedpe.strand_b):
            return False

        if (a_bedpe.chrom_a != b_bedpe.chrom_a
            or a_bedpe.start_a - max_distance > b_bedpe.end_a
            or a_bedpe.end_a + max_distance < b_bedpe.start_a):
            return False

        if (a_bedpe.chrom_b != b_bedpe.chrom_b
            or a_bedpe.start_b - max_distance > b_bedpe.end_b
            or a_bedpe.end_b + max_distance < b_bedpe.start_b):
            return False
        else:
            a_bedpe.cohort_vars[b_bedpe.id]=b_bedpe.af
            return True
    else:
        return False
def varLookup(aFile, bFile,bedpe_out, max_distance,pass_prefix,cohort_name):
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
        bList.append(Bedpe(bLine))
    
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
            a = Bedpe(aLine)
            for b in bList:
                add(a,b,max_distance)
            bedpe_out.write(a.get_var_string(cohort_name))

def description():
    return 'Look for variants common between two bedpe files'

def add_arguments_to_parser(parser):
    parser.add_argument('-d', '--distance', type=int, dest='max_distance', default=50, help='max separation distance (bp) of adjacent loci between bedpe files [50]')
    parser.add_argument("-a", "--aFile", dest="aFile", help="Pruned merged bedpe (A file) or standard input (-a stdin).")
    parser.add_argument("-b", "--bFile", dest="bFile", help="Pruned merged bedpe (B file) (-b stdin). For pruning use svtools prune")
    parser.add_argument("-c", "--cohort", dest='cohort_name', default=None, help="Cohort name to add information of matching variants (default:bFile)")                    
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help='Output BEDPE to write (default: stdout)')
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
    	varLookup(args.aFile, args.bFile,args.output, args.max_distance,pass_prefix,args.cohort_name)
    except IOError as err:
    	sys.stderr.write("IOError " + str(err) + "\n");

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
