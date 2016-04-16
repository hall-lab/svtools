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
        self.coordinate_buffer = 500
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

    def distant_coordinates(self, strand, start, stop):
        '''
        Transform coordinates if SV is too large
        '''
        if strand == '+':
            stop += self.coordinate_buffer
        else:
            start -= self.coordinate_buffer
        return (start, stop)

    def distant_block_sizes(self, strand, start, stop):
        '''
        Return tuple of chunk sizes
        '''
        rv = (stop - start, 1)
        if strand == '+':
            return rv
        else:
            return rv[::-1] #reverse

    def distant_block_starts(self, strand, start, stop):
        '''
        Return tuple of chunk sizes
        '''
        second_start = self.coordinate_buffer
        if strand == '+':
            second_start += (stop - start) - 1
        return (0, second_start)

    @staticmethod
    def create_line(chrom, start, end, name, score, strand, color, size_tuple=None, start_tuple=None):
        '''
        Create a blockedBed line
        '''
        fields = [
            chrom, 
            start, 
            end, 
            name, 
            score, 
            strand, 
            start, 
            end, 
            color
            ]
        if size_tuple is not None and start_tuple is not None:
            fields += [
                    '2',
                    ','.join(map(str, size_tuple)),
                    ','.join(map(str,start_tuple))
            ]
        return '\t'.join(map(str,fields))

    def convert(self, bedpe):
        '''
        Convert Bedpe line to a Bed12 line
        '''
        span = abs(bedpe.e2 - bedpe.s1)
        color = self.get_color(bedpe.svtype, span)

        output_lines = list()
        if (bedpe.svtype != 'BND'):
            name = self.bed12_name(bedpe.svtype, bedpe.name, bedpe.af)
            if span <= self.max_dist:
                output_lines.append(self.create_line(
                    bedpe.c1, 
                    bedpe.s1, 
                    bedpe.e2, 
                    name, 
                    bedpe.score, 
                    '+',
                    color, 
                    (bedpe.e1 - bedpe.s1, bedpe.e2 - bedpe.s2),
                    (0, bedpe.s2 - bedpe.s1)))
            else:
                s1, e1 = self.distant_coordinates(bedpe.o1, bedpe.s1, bedpe.e1)
                size_tuple1 = self.distant_block_sizes(bedpe.o1, bedpe.s1, bedpe.e1)
                start_tuple1 = self.distant_block_starts(bedpe.o1, bedpe.s1, bedpe.e1)
                output_lines.append(self.create_line(
                    bedpe.c1, 
                    s1, 
                    e1, 
                    name, 
                    bedpe.score, 
                    bedpe.o1, 
                    color, 
                    size_tuple1,
                    start_tuple1))
                s2, e2 = self.distant_coordinates(bedpe.o2, bedpe.s2, bedpe.e2)
                size_tuple2 = self.distant_block_sizes(bedpe.o2, bedpe.s2, bedpe.e2)
                start_tuple2 = self.distant_block_starts(bedpe.o2, bedpe.s2, bedpe.e2)
                output_lines.append(self.create_line(
                    bedpe.c2,
                    s2,
                    e2,
                    name,
                    bedpe.score,
                    bedpe.o2,
                    color,
                    size_tuple2,
                    start_tuple2))
        else:
            name = self.bed12_name(bedpe.svtype, bedpe.name, bedpe.af, (bedpe.o1, bedpe.o2))
            output_lines.append(self.create_line(
                bedpe.c1,
                bedpe.s1,
                bedpe.e1,
                name, 
                bedpe.score,
                bedpe.o1,
                color))
            output_lines.append(self.create_line(
                bedpe.c2,
                bedpe.s2,
                bedpe.e2,
                name, 
                bedpe.score,
                bedpe.o2,
                color))
        return output_lines

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
