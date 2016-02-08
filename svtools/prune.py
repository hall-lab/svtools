import argparse, sys , re
from argparse import RawTextHelpFormatter
from collections import Counter
from svtools.bedpe import Bedpe
import svtools.utils as su

class Cluster(object):
    def __init__(self):
        self.id = None
        self.elements = [None]
        self.chrom_a = None
        self.min_a = float("inf")
        self.max_a = 0
        self.chrom_b = None
        self.min_b = float("inf")
        self.max_b = 0
        self.size = 0
        self.strand_a = ''
        self.strand_b = ''
        self.sv_event = ''
        self.filter = 0

    # check whether a bedpe object is addable to this
    # cluster given the max_distance
    def can_add(self, bedpe, max_distance):
        if self.size == 0:
           return True
        if self.sv_event ==  bedpe.svtype:
            if (self.strand_a != bedpe.o1
                or self.strand_b != bedpe.o2):
                return False

            if (self.chrom_a != bedpe.c1
                or self.min_a - max_distance > bedpe.e1
                or self.max_a + max_distance < bedpe.s1):
                return False

            if (self.chrom_b != bedpe.c2
                or self.min_b - max_distance > bedpe.e2
                or self.max_b + max_distance < bedpe.s2):
                return False
            else:
                return True
        else:
            return False

    def add(self, bedpe, max_distance, eval_param):
        if eval_param == None or eval_param.lower() == 'af':
            if  bedpe.af != '.' and bedpe.af > self.filter:
                #First node represents best variant
                self.elements[0] = bedpe  
        self.size += 1
        self.sv_event=bedpe.svtype
        self.filter = max(self.filter,bedpe.af) 
        self.chrom_a = bedpe.c1
        self.min_a = min(self.min_a, bedpe.s1)
        self.max_a = max(self.max_a, bedpe.e1)
        self.chrom_b = bedpe.c2
        self.min_b = min(self.min_b, bedpe.s2)
        self.max_b = max(self.max_b, bedpe.e2)
        self.strand_a = bedpe.o1
        self.strand_b = bedpe.o2
                 
    def get_cluster_string(self):
        # FIXME Should move BEDPE -> string stuff to bedpe class
        b = self.elements[0]
        return '\t'.join([
            b.c1,
            str(b.s1),
            str(b.e1),
            b.c2,
            str(b.s2),
            str(b.e2),
            b.name,
            str(b.score),
            b.o1,
            b.o2,
            b.svtype,
            b.filter,
            b.misc[0],
            b.info2
            ] + b.misc[1:])


# prints and removes clusters from cluster_list that are beyond
# distance window
def prune(cluster_list, bedpe, max_distance, min_cluster_size,print_ineligible,bedpe_out):
    new_cluster_list = []
    global cluster_counter
    for cluster in cluster_list:
        # cluster is beyond updatable window:
        if (bedpe is None   
            or cluster.chrom_a != bedpe.c1
            or cluster.min_a - max_distance > bedpe.e1
            or cluster.max_a + max_distance < bedpe.s1):
            
            # print the cluster if eligible
            if (cluster.size >= min_cluster_size or print_ineligible):
                bedpe_out.write(cluster.get_cluster_string() + '\n')
                cluster_counter += 1
                cluster.id = cluster_counter
            else:
                new_cluster_list.append(cluster)
        # cluster is still within updatable window,
        # leave it in the cluster list
        else:
            new_cluster_list.append(cluster)

    return new_cluster_list

# primary function
def cluster_bedpe(in_file,max_distance,eval_param,bedpe_out,is_sorted):
    
    global min_cluster_size,cluster_counter
    #minimum cluster size to report feature
    min_cluster_size = 1
    # Inititalize the number of clusters
    cluster_counter = 0
    # line number
    line_counter = 0
    # a list of clusters in the buffer
    cluster_list = []
    #Flag to print features which are not 
    #eligible according to min_cluster_size(<2)
    print_ineligible = 0    
    #Print headers at top
    in_header = True
    for line in in_file:
        if line.startswith('#') and in_header:
            bedpe_out.write(line)
            continue
        in_header=False
        line_counter += 1
        bedpe = Bedpe(line.rstrip().split('\t'))
        if bedpe.af is None:
            sys.stderr.write('No allele frequency for variant found. This tool requires allele frequency information to function. Please add with svtools afreq and rerun\n')
            sys.exit(1)

        matched_clusters=[]
        for cluster in cluster_list:
            if cluster.can_add(bedpe, max_distance):
                cluster.add(bedpe, max_distance, eval_param)
                matched_clusters.append(cluster)
        if len(matched_clusters) == 0:
            new_cluster = Cluster()
            new_cluster.add(bedpe, max_distance, eval_param)
            cluster_list.append(new_cluster)
        else:
            if len(matched_clusters) > 1:
                i = 0
                matched_cluster_pruned = False
                while i < (len(matched_clusters)-1):
                    j=i+1
                    while j < len(matched_clusters):
                        if matched_clusters[i].can_add(matched_clusters[j].elements[0],max_distance):
                            matched_clusters[i].add(matched_clusters[j].elements[0], max_distance, eval_param)
                            matched_cluster_pruned = True
                            del matched_clusters[j]
                        j+=1
                    i+=1        
                if matched_cluster_pruned == True:
                        cluster_list = [cluster for cluster in cluster_list if cluster not in matched_clusters]
        #prune and print eligible clusters
        if line_counter % 1000 == 0 and is_sorted:
            cluster_list = prune(cluster_list,
                                 bedpe,
                                 max_distance,
                                 min_cluster_size,
                                 print_ineligible,
                                 bedpe_out)

    # prune the cluster and print eligible
    # features
    ##
    print_ineligible = True 
    cluster_list = prune(cluster_list,
                         None,
                         max_distance,
                         min_cluster_size,
                         print_ineligible,
                         bedpe_out)
        
    # close the input file
    in_file.close()
    return

def description():
    return 'Cluster and prune a BEDPE file by position based on allele frequency'

def add_arguments_to_parser(parser):
    parser.add_argument('-d', '--distance', type=int, dest='max_distance', default=50, help='max separation distance (bp) of adjacent loci in cluster [50]')
    parser.add_argument('-e', '--eval-param', help='evaluating parameter for choosing best bedpe in a cluster(e.g. af=AlleleFrequency default:af)')
    parser.add_argument('-s', '--is-sorted', action='store_true', help='specify if an input file is sorted. Sort with svtools bedpesort. (default=False)')
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=None, help='BEDPE file to read. If \'-\' or absent then defaults to stdin.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help='Output bedpe to write (default: stdout)')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.input) as stream:
        cluster_bedpe(stream,
                  args.max_distance,
                  args.eval_param,
                  args.output,
                  args.is_sorted)

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
