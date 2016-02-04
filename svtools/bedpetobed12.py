import sys
import argparse
from svtools.bedpe import Bedpe

def bedpeToBlockedBed(bedpe, dist):

    if (bedpe.svtype == "DEL") and (abs(bedpe.e2-bedpe.s1) <= dist): color = "153,0,0"    # deletion breakpoints are red
    elif (bedpe.svtype == "DUP") and (abs(bedpe.e2-bedpe.s1) <= dist) : color = "0,102,0"  # duplication breakpoints are green
    elif (bedpe.svtype == "INV") and (abs(bedpe.e2-bedpe.s1) <= dist) : color = "0,51,204"  # inversion breakpoints are blue
    elif (bedpe.svtype == "BND") or (abs(bedpe.e2-bedpe.s1) > dist): color = "204,204,204" # distant breakpoints are gray
    elif abs(bedpe.e2-bedpe.s1) <= dist:color="128,0,128" 
    if (bedpe.svtype != "BND") and (abs(bedpe.e2-bedpe.s1) <= dist):
        if bedpe.af is not None:
            print '\t'.join(map(str, [bedpe.c1,bedpe.s1,bedpe.e2,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name,';AF=',bedpe.af])),bedpe.score,'+',\
                            bedpe.s1,bedpe.e2,color,'2',','.join(map(str,[bedpe.e1-bedpe.s1,bedpe.e2-bedpe.s2])), ','.join(map(str,['0', bedpe.s2 - bedpe.s1]))]))
        else:
            print '\t'.join(map(str, [bedpe.c1,bedpe.s1,bedpe.e2,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name])),bedpe.score,'+',\
                            bedpe.s1,bedpe.e2,color,'2',','.join(map(str,[bedpe.e1-bedpe.s1,bedpe.e2-bedpe.s2])), ','.join(map(str,['0', bedpe.s2 - bedpe.s1]))]))
    # intrachromosomals that exceed dist
    elif (bedpe.svtype != "BND") and (abs(bedpe.e2-bedpe.s1) > dist):
        if bedpe.o1 == "+":
            if bedpe.af is not None:
                print '\t'.join(map(str, [bedpe.c1,bedpe.s1,bedpe.e1+500,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name,';AF=',bedpe.af])),bedpe.score,'+',\
                                bedpe.s1,bedpe.e1+500,color,'2',','.join(map(str,[bedpe.e1-bedpe.s1,1])), ','.join(map(str,[0, bedpe.e1 - bedpe.s1+500]))]))
            else:
                print '\t'.join(map(str, [bedpe.c1,bedpe.s1,bedpe.e1+500,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name])),bedpe.score,'+',\
                                bedpe.s1,bedpe.e1+500,color,'2',','.join(map(str,[bedpe.e1-bedpe.s1,1])), ','.join(map(str,[0, bedpe.e1 - bedpe.s1+500]))]))
        if bedpe.o1 == "-":
            if bedpe.af is not None:
                print '\t'.join(map(str, [bedpe.c1,bedpe.s1-500,bedpe.e1,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name,';AF=',bedpe.af])),bedpe.score,'-',\
                                bedpe.s1-500,bedpe.e1,color,'2',','.join(map(str,[1,bedpe.e1-bedpe.s1])), ','.join(map(str,[0, 500]))]))
            else:
                print '\t'.join(map(str, [bedpe.c1,bedpe.s1-500,bedpe.e1,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name])),bedpe.score,'-',\
                                bedpe.s1-500,bedpe.e1,color,'2',','.join(map(str,[1,bedpe.e1-bedpe.s1])), ','.join(map(str,[0, 500]))]))
        if bedpe.o2 == "+":
            if bedpe.af is not None:
                print '\t'.join(map(str, [bedpe.c2,bedpe.s2,bedpe.e2+500,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name,';AF=',bedpe.af])),bedpe.score,'+',\
                                bedpe.s2,bedpe.e2+500,color,'2',','.join(map(str,[bedpe.e2-bedpe.s2,1])), ','.join(map(str,[0, bedpe.e2-bedpe.s2+499]))]))
            else:
                print '\t'.join(map(str, [bedpe.c2,bedpe.s2,bedpe.e2+500,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name])),bedpe.score,'+',\
                                bedpe.s2,bedpe.e2+500,color,'2',','.join(map(str,[bedpe.e2-bedpe.s2,1])), ','.join(map(str,[0, bedpe.e2-bedpe.s2+499]))]))
            
        if bedpe.o2 == "-":
            print bedpe.c2 + "\t" + str(bedpe.s2-500) + "\t" + str(bedpe.e2) + "\t"  + bedpe.svtype + "_" + bedpe.name + \
            "\t" + str(bedpe.score) + "\t" + \
            "-" + "\t" + str(bedpe.s2-500) + "\t" + str(bedpe.e2) + "\t" + color + "\t" + "2" + "\t" + \
            "1" + "," + str(bedpe.e2 - bedpe.s2) + "\t" + \
            "0," + str(500)
            if bedpe.af is not None:
                print '\t'.join(map(str, [bedpe.c2,bedpe.s2-500,bedpe.e2,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name,';AF=',bedpe.af])),bedpe.score,'-',\
                                bedpe.s2-500,bedpe.e2,color,'2',','.join(map(str,[1,bedpe.e2-bedpe.s2,1])), ','.join(map(str,[0,500]))]))
            else:
                print '\t'.join(map(str, [bedpe.c2,bedpe.s2-500,bedpe.e2,''.join(map(str,[bedpe.svtype,';ID=',bedpe.name])),bedpe.score,'-',\
                                bedpe.s2-500,bedpe.e2,color,'2',','.join(map(str,[1,bedpe.e2-bedpe.s2,1])), ','.join(map(str,[0,500]))]))
    # BNDS:	
    elif (bedpe.svtype == "BND"):
        if bedpe.af is not None:
            print '\t'.join([bedpe.c1,str(bedpe.s1),str(bedpe.e1),''.join(["BND;ID=",bedpe.name,";AF=",bedpe.af,";STR=",bedpe.o1,bedpe.o2]), \
                            str(bedpe.score),bedpe.o1,str(bedpe.s1),str(bedpe.e1),"204,204,204"]) 
            print '\t'.join([bedpe.c2,str(bedpe.s2),str(bedpe.e2),''.join(["BND;ID=",bedpe.name,";AF=",bedpe.af,";STR=",bedpe.o1,bedpe.o2]),\
                            str(bedpe.score),bedpe.o2,str(bedpe.s2),str(bedpe.e2),"204,204,204"])
        else:
            print '\t'.join([bedpe.c1,str(bedpe.s1),str(bedpe.e1),''.join(["BND;ID=",bedpe.name,";STR=",bedpe.o1,bedpe.o2]),\
                            str(bedpe.score),bedpe.o1,str(bedpe.s1),str(bedpe.e1),"204,204,204"])
            print '\t'.join([bedpe.c2,str(bedpe.s2),str(bedpe.e2),''.join(["BND;ID=",bedpe.name,";STR=",bedpe.o1,bedpe.o2]),\
                            str(bedpe.score),bedpe.o2,str(bedpe.s2),str(bedpe.e2),"204,204,204"]) 
            
def processBEDPE(bedpeFile, name, dist):
    #Process the BEDPE file and convert each entry to SAM.
    if name is not None or bedpeFile == "stdin":
        writeTrackName(name)
    elif bedpeFile != "stdin":
        writeTrackName(bedpeFile.name)    
    if bedpeFile == "stdin":		
        for line in sys.stdin:
            # ignore header
            if line[0] == "#":
                continue
            lineList = line.strip().split()
            if (len(lineList) > 0):
                bedpe = Bedpe(lineList)
                bedpeToBlockedBed(bedpe, dist)
    else:
         for line in open(bedpeFile, 'r'):
             # ignore header
            if line[0] == "#":
                 continue
            lineList = line.strip().split()
            if (len(lineList) > 0):
                bedpe = Bedpe(lineList)
                bedpeToBlockedBed(bedpe, dist)



def writeTrackName(name):
    print "track name=" + name + " itemRgb=On"

def description():
    return 'converts BEDPE to BED12 format for viewing in IGV or the UCSC browser'

def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--bedpe', help='BEDPE input file', required=True)
    parser.add_argument('-n', '--name', default='BEDPE', help="The name of the track. Default is 'BEDPE'")
    parser.add_argument('-d', '--maxdist', dest='dist', default=1000000, type=int, help='The minimum distance for drawing intrachromosomal features as if they are interchromosomal (i.e., without a line spanning the two footprints). Default is 1Mb.')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    processBEDPE(args.bedpe, args.name, args.dist)

if __name__ == "__main__":
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
