import sys
import argparse
from svtools.bedpe import Bedpe
import svtools.utils as su

class BedpetoBlockedBedConverter(object):
    '''
    Class to convert Bedpe lines to BlockedBed (Bed12)
    '''
    def __init__(self, name, max_dist):
        '''
        Initialize parameters
        '''
        self.max_dist = max_dist
        self.name = name
        self.distant_color = '204,204,204' #gray
        self.unknown_close_color = '128,0,128'
        self.color_table = { 
                'DEL': '153,0,0', #red
                'DUP': '0,102,0', #green
                'INV': '0,51,204', #blue
                'BND': self.distant_color
                }
    
    def track_name(self):
        '''
        Return the track name for the output bedpe. Assumes name is a valid string
        '''
        return "track name=" + self.name + " itemRgb=On\n"

    def get_color(self, svtype, span):
        '''
        Get the color for the Bed12 entries
        '''
        if span <= self.max_dist:
            return self.color_table.get(svtype, self.unknown_close_color)
        else:
            return self.distant_color

    @staticmethod
    def bed12_name(svtype, bedpe_name, allele_frequency, strands=None):
        '''
        Construct the name for the Bed12 entry.
        '''
        base_string = '{0};ID={1}'.format(svtype, bedpe_name)
        if allele_frequency is not None:
            base_string += ';AF={0}'.format(allele_frequency)
        if strands is not None:
            base_string += ';STR={0}'.format(''.join(strands))
        return base_string

    def convert(self, bedpe):
        '''
        Convert Bedpe line to a Bed12 line
        '''
        span = abs(bedpe.e2 - bedpe.s1)
        color = self.get_color(bedpe.svtype, span)

        output_lines = list()
        if (bedpe.svtype != "BND") and (span <= self.max_dist):
            name = self.bed12_name(bedpe.svtype, bedpe.name, bedpe.af)
            output_lines.append(
                    '\t'.join(
                        map(
                            str, 
                            [
                                bedpe.c1,
                                bedpe.s1,
                                bedpe.e2,
                                name,
                                bedpe.score,
                                '+',
                                bedpe.s1,
                                bedpe.e2,
                                color,
                                '2',
                                ','.join(map(str,[bedpe.e1-bedpe.s1,bedpe.e2-bedpe.s2])), 
                                ','.join(map(str,['0', bedpe.s2 - bedpe.s1]))
                                ]
                            )
                        )
                    )
        # intrachromosomals that exceed dist
        elif (bedpe.svtype != "BND") and  (span > self.max_dist):
            name = self.bed12_name(bedpe.svtype, bedpe.name, bedpe.af)
            if bedpe.o1 == "+":
                output_lines.append(
                        '\t'.join(
                            map(
                                str,
                                [
                                    bedpe.c1,
                                    bedpe.s1,
                                    bedpe.e1+500,
                                    name,
                                    bedpe.score,
                                    '+',
                                    bedpe.s1,
                                    bedpe.e1+500,
                                    color,
                                    '2',
                                    ','.join(map(str,[bedpe.e1-bedpe.s1,1])),
                                    ','.join(map(str,[0, bedpe.e1 - bedpe.s1+500]))
                                    ]
                                )
                            )
                        )
            if bedpe.o1 == "-":
                output_lines.append(
                        '\t'.join(
                            map(
                                str,
                                [
                                    bedpe.c1,
                                    bedpe.s1-500,
                                    bedpe.e1,
                                    name,
                                    bedpe.score,
                                    '-',
                                    bedpe.s1-500,
                                    bedpe.e1,
                                    color,
                                    '2',
                                    ','.join(map(str,[1,bedpe.e1-bedpe.s1])),
                                    ','.join(map(str,[0, 500]))]
                                )
                            )
                        )
            if bedpe.o2 == "+":
                output_lines.append(
                        '\t'.join(
                            map(
                                str, 
                                [
                                    bedpe.c2,
                                    bedpe.s2,
                                    bedpe.e2+500,
                                    name,
                                    bedpe.score,
                                    '+',
                                    bedpe.s2,
                                    bedpe.e2+500,
                                    color,
                                    '2',
                                    ','.join(map(str,[bedpe.e2-bedpe.s2,1])), 
                                    ','.join(map(str,[0, bedpe.e2-bedpe.s2+499]))])))
            if bedpe.o2 == "-":
                output_lines.append(
                        '\t'.join(
                            map(
                                str,
                                [
                                    bedpe.c2,
                                    bedpe.s2-500,
                                    bedpe.e2,
                                    name,
                                    bedpe.score,
                                    '-',
                                    bedpe.s2-500,
                                    bedpe.e2,
                                    color,
                                    '2',
                                    ','.join(map(str,[1,bedpe.e2-bedpe.s2])),
                                    ','.join(map(str,[0,500]))
                                    ]
                                )
                            )
                        )

        # BNDS:	
        elif (bedpe.svtype == "BND"):
            name = self.bed12_name(bedpe.svtype, bedpe.name, bedpe.af, (bedpe.o1, bedpe.o2))
            output_lines.append(
                    '\t'.join(
                        map(
                            str,
                            [
                                bedpe.c1,
                                bedpe.s1,
                                bedpe.e1,
                                name, 
                                bedpe.score,
                                bedpe.o1,
                                bedpe.s1,
                                bedpe.e1,
                                color
                                ]
                            )
                        )
                    )
            output_lines.append(
                    '\t'.join(
                        map(
                            str,
                            [
                                bedpe.c2,
                                bedpe.s2,
                                bedpe.e2,
                                name, 
                                bedpe.score,
                                bedpe.o2,
                                bedpe.s2,
                                bedpe.e2,
                                color
                                ]
                            )
                        )
                    )
        return output_lines

def bedpeToBlockedBed(bedpe, dist, output_handle=sys.stdout):

    if (bedpe.svtype == "DEL") and (abs(bedpe.e2-bedpe.s1) <= dist): color = "153,0,0"    # deletion breakpoints are red
    elif (bedpe.svtype == "DUP") and (abs(bedpe.e2-bedpe.s1) <= dist) : color = "0,102,0"  # duplication breakpoints are green
    elif (bedpe.svtype == "INV") and (abs(bedpe.e2-bedpe.s1) <= dist) : color = "0,51,204"  # inversion breakpoints are blue
    elif (bedpe.svtype == "BND") or (abs(bedpe.e2-bedpe.s1) > dist): color = "204,204,204" # distant breakpoints are gray
    elif abs(bedpe.e2-bedpe.s1) <= dist:color="128,0,128" 
    output_lines = list()
    if (bedpe.svtype != "BND") and (abs(bedpe.e2-bedpe.s1) <= dist):
        if bedpe.af is not None:
            output_lines.append('\t'.join(map(str, [bedpe.c1,bedpe.s1,bedpe.e2,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name,';AF=',bedpe.af])),bedpe.score,'+',\
                            bedpe.s1,bedpe.e2,color,'2',','.join(map(str,[bedpe.e1-bedpe.s1,bedpe.e2-bedpe.s2])), ','.join(map(str,['0', bedpe.s2 - bedpe.s1]))])))
        else:
            output_lines.append('\t'.join(map(str, [bedpe.c1,bedpe.s1,bedpe.e2,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name])),bedpe.score,'+',\
                            bedpe.s1,bedpe.e2,color,'2',','.join(map(str,[bedpe.e1-bedpe.s1,bedpe.e2-bedpe.s2])), ','.join(map(str,['0', bedpe.s2 - bedpe.s1]))])))
    # intrachromosomals that exceed dist
    elif (bedpe.svtype != "BND") and (abs(bedpe.e2-bedpe.s1) > dist):
        if bedpe.o1 == "+":
            if bedpe.af is not None:
                output_lines.append('\t'.join(map(str, [bedpe.c1,bedpe.s1,bedpe.e1+500,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name,';AF=',bedpe.af])),bedpe.score,'+',\
                                bedpe.s1,bedpe.e1+500,color,'2',','.join(map(str,[bedpe.e1-bedpe.s1,1])), ','.join(map(str,[0, bedpe.e1 - bedpe.s1+500]))])))
            else:
                output_lines.append('\t'.join(map(str, [bedpe.c1,bedpe.s1,bedpe.e1+500,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name])),bedpe.score,'+',\
                                bedpe.s1,bedpe.e1+500,color,'2',','.join(map(str,[bedpe.e1-bedpe.s1,1])), ','.join(map(str,[0, bedpe.e1 - bedpe.s1+500]))])))
        if bedpe.o1 == "-":
            if bedpe.af is not None:
                output_lines.append('\t'.join(map(str, [bedpe.c1,bedpe.s1-500,bedpe.e1,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name,';AF=',bedpe.af])),bedpe.score,'-',\
                                bedpe.s1-500,bedpe.e1,color,'2',','.join(map(str,[1,bedpe.e1-bedpe.s1])), ','.join(map(str,[0, 500]))])))
            else:
                output_lines.append('\t'.join(map(str, [bedpe.c1,bedpe.s1-500,bedpe.e1,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name])),bedpe.score,'-',\
                                bedpe.s1-500,bedpe.e1,color,'2',','.join(map(str,[1,bedpe.e1-bedpe.s1])), ','.join(map(str,[0, 500]))])))
        if bedpe.o2 == "+":
            if bedpe.af is not None:
                output_lines.append('\t'.join(map(str, [bedpe.c2,bedpe.s2,bedpe.e2+500,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name,';AF=',bedpe.af])),bedpe.score,'+',\
                                bedpe.s2,bedpe.e2+500,color,'2',','.join(map(str,[bedpe.e2-bedpe.s2,1])), ','.join(map(str,[0, bedpe.e2-bedpe.s2+499]))])))
            else:
                output_lines.append('\t'.join(map(str, [bedpe.c2,bedpe.s2,bedpe.e2+500,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name])),bedpe.score,'+',\
                                bedpe.s2,bedpe.e2+500,color,'2',','.join(map(str,[bedpe.e2-bedpe.s2,1])), ','.join(map(str,[0, bedpe.e2-bedpe.s2+499]))])))
            
        if bedpe.o2 == "-":
            if bedpe.af is not None:
                output_lines.append('\t'.join(map(str, [bedpe.c2,bedpe.s2-500,bedpe.e2,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name,';AF=',bedpe.af])),bedpe.score,'-',\
                                bedpe.s2-500,bedpe.e2,color,'2',','.join(map(str,[1,bedpe.e2-bedpe.s2])), ','.join(map(str,[0,500]))])))
            else:
                output_lines.append('\t'.join(map(str, [bedpe.c2,bedpe.s2-500,bedpe.e2,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name])),bedpe.score,'-',\
                                bedpe.s2-500,bedpe.e2,color,'2',','.join(map(str,[1,bedpe.e2-bedpe.s2])), ','.join(map(str,[0,500]))])))
    # BNDS:	
    elif (bedpe.svtype == "BND"):
        if bedpe.af is not None:
            output_lines.append('\t'.join([bedpe.c1,str(bedpe.s1),str(bedpe.e1),''.join(["BND;ID=",bedpe.name,";AF=",bedpe.af,";STR=",bedpe.o1,bedpe.o2]), \
                            str(bedpe.score),bedpe.o1,str(bedpe.s1),str(bedpe.e1),"204,204,204"]))
            output_lines.append('\t'.join([bedpe.c2,str(bedpe.s2),str(bedpe.e2),''.join(["BND;ID=",bedpe.name,";AF=",bedpe.af,";STR=",bedpe.o1,bedpe.o2]),\
                            str(bedpe.score),bedpe.o2,str(bedpe.s2),str(bedpe.e2),"204,204,204"]))
        else:
            output_lines.append('\t'.join([bedpe.c1,str(bedpe.s1),str(bedpe.e1),''.join(["BND;ID=",bedpe.name,";STR=",bedpe.o1,bedpe.o2]),\
                            str(bedpe.score),bedpe.o1,str(bedpe.s1),str(bedpe.e1),"204,204,204"]))
            output_lines.append('\t'.join([bedpe.c2,str(bedpe.s2),str(bedpe.e2),''.join(["BND;ID=",bedpe.name,";STR=",bedpe.o1,bedpe.o2]),\
                            str(bedpe.score),bedpe.o2,str(bedpe.s2),str(bedpe.e2),"204,204,204"]))
    if len(output_lines) > 0:
        output_handle.write('\n'.join(output_lines) + '\n')
            
def processBEDPE(bedpe_stream, name, dist, output_handle=sys.stdout):
    #Process the BEDPE file and convert each entry to SAM.
    converter = BedpetoBlockedBedConverter(name, dist)
    output_handle.write(converter.track_name())
    for line in bedpe_stream:
        # ignore header
        if line[0] == "#":
            continue
        lineList = line.rstrip().split('\t')
        if lineList:
            bedpe = Bedpe(lineList)
            output_handle.write('\n'.join(converter.convert(bedpe)) + '\n')

def writeTrackName(name, output_handle=sys.stdout):
    output_handle.write("track name=" + name + " itemRgb=On\n")

def description():
    return 'convert a BEDPE file to BED12 format for viewing in IGV or the UCSC browser'

def epilog():
    return 'The input BEDPE file may be gzipped. If the input file is omitted then input is read from stdin. Output is written to stdout.' 

def add_arguments_to_parser(parser):
    parser.add_argument('-b', '--bedpe', metavar='<BEDPE>', default=None, help='BEDPE input file')
    parser.add_argument('-n', '--name', metavar='<STRING>', default='BEDPE', help="The name of the track. Default is 'BEDPE'")
    parser.add_argument('-d', '--maxdist', metavar='<INT>', dest='dist', default=1000000, type=int, help='The minimum distance for drawing intrachromosomal features as if they are interchromosomal (i.e., without a line spanning the two footprints). Default is 1Mb.')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description(), epilog=epilog())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.bedpe) as stream:
        processBEDPE(stream, args.name, args.dist)

if __name__ == "__main__":
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
