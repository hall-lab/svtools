import sys
import argparse
from svtools.bedpe import Bedpe
import svtools.utils as su

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
            output_lines.append(bedpe.c2 + "\t" + str(bedpe.s2-500) + "\t" + str(bedpe.e2) + "\t"  + bedpe.svtype + "_" + bedpe.name + \
            "\t" + str(bedpe.score) + "\t" + \
            "-" + "\t" + str(bedpe.s2-500) + "\t" + str(bedpe.e2) + "\t" + color + "\t" + "2" + "\t" + \
            "1" + "," + str(bedpe.e2 - bedpe.s2) + "\t" + \
            "0," + str(500))
            if bedpe.af is not None:
                output_lines.append('\t'.join(map(str, [bedpe.c2,bedpe.s2-500,bedpe.e2,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name,';AF=',bedpe.af])),bedpe.score,'-',\
                                bedpe.s2-500,bedpe.e2,color,'2',','.join(map(str,[1,bedpe.e2-bedpe.s2,1])), ','.join(map(str,[0,500]))])))
            else:
                output_lines.append('\t'.join(map(str, [bedpe.c2,bedpe.s2-500,bedpe.e2,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name])),bedpe.score,'-',\
                                bedpe.s2-500,bedpe.e2,color,'2',','.join(map(str,[1,bedpe.e2-bedpe.s2,1])), ','.join(map(str,[0,500]))])))
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
    if name is not None or bedpe_stream == sys.stdin:
        writeTrackName(name, output_handle)
    else:
        writeTrackName(bedpe_stream.name, output_handle)    
    for line in bedpe_stream:
        # ignore header
        if line[0] == "#":
            continue
        lineList = line.strip().split()
        if (len(lineList) > 0):
                bedpe = Bedpe(lineList)
                bedpeToBlockedBed(bedpe, dist, output_handle)

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
