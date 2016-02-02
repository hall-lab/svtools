#!/usr/bin/env python

import argparse, sys , re
from argparse import RawTextHelpFormatter
from collections import Counter

__author__ = "Abhijit Badve"
__version__ = "$Revision: 0.0.1 $"
___date__ = "$Date: 2015-08-25 12:54 $"

# --------------------------------------
# define functions
def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
clusterBedpe.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: cluster a BEDPE file by position based on their allele frequency")
    parser.add_argument('-d', '--distance',
                        metavar='INT', type=int,
                        dest='max_distance',
                        required=False,
                        default=50,
                        help='max separation distance (bp) of adjacent loci in cluster [50]')
    parser.add_argument('-e', '--eval_param',
                        metavar='string', type=str,
                        required=False,
                        dest='eval_param',
                        help='evaluating parameter for choosing best bedpe in a cluster(e.g. af=AlleleFrequency default:af)')
    parser.add_argument('-s', '--is_sorted',
                        required=False,
                        action='store_true',
                        help='specifying if an input file is sorted (default=False)\n(use following command to sort: \'cat FILE | sort -k1,1V -k2,2n -k3,3n -k4,4V -k5,5n -k6,6\')') 
                                            
    parser.add_argument('input', nargs='?',
                        type=argparse.FileType('r'),
                        default=None,
                        help='BEDPE file to read. If \'-\' or absent then defaults to stdin.')
    parser.add_argument('-o', '--output', 
                        type=argparse.FileType('w'),
                        default=sys.stdout, 
                        help='Output bedpe to write (default: stdout)')
    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input = sys.stdin

    # send back the user input
    return args

class Bedpe(object):
    def __init__(self, line):
        v = line.rstrip().split("\t")
        self.chrom_a = v[0]
        self.start_a = int(v[1])
        self.end_a = int(v[2])
        self.chrom_b = v[3]
        self.start_b = int(v[4])
        self.end_b = int(v[5])
        self.id = v[6]
        self.af=filter(lambda x: x if x.startswith('AF=') else None, v[12].split(";"))
        if self.af is not None and len(self.af) == 1:
           self.af=''.join(self.af).replace('AF=','')
        else: 
           print "No allele frequency for variants found. Input a valid file"
           sys.exit(0)
        
        self.sv_event=filter(lambda x: 'SVTYPE=' in x, v[12].split(";"))
        if len(self.sv_event) == 1:
            self.sv_event=''.join(self.sv_event).replace('SVTYPE=','')
        else:
            print "No SVTYPE Field found. Input a valid file"
            sys.exit(0)
        # self.sv_event=re.search('SVTYPE=(.+?);',v[12]).group(1);
        self.vec = '\t'.join(v)
        try:
            self.strand_a = v[8]
            self.strand_b = v[9]
        except IndexError:
            self.strand_a = ''
            self.strand_b = ''
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
        if self.sv_event ==  bedpe.sv_event:
            if (self.strand_a != bedpe.strand_a
                or self.strand_b != bedpe.strand_b):
                return False

            if (self.chrom_a != bedpe.chrom_a
                or self.min_a - max_distance > bedpe.end_a
                or self.max_a + max_distance < bedpe.start_a):
                return False

            if (self.chrom_b != bedpe.chrom_b
                or self.min_b - max_distance > bedpe.end_b
                or self.max_b + max_distance < bedpe.start_b):
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
        self.sv_event=bedpe.sv_event 
        self.filter = max(self.filter,bedpe.af) 
        self.chrom_a = bedpe.chrom_a
        self.min_a = min(self.min_a, bedpe.start_a)
        self.max_a = max(self.max_a, bedpe.end_a)
        self.chrom_b = bedpe.chrom_b
        self.min_b = min(self.min_b, bedpe.start_b)
        self.max_b = max(self.max_b, bedpe.end_b)
        self.strand_a = bedpe.strand_a
        self.strand_b = bedpe.strand_b
                 
    def get_cluster_string(self):
           return self.elements[0].vec


# prints and removes clusters from cluster_list that are beyond
# distance window
def prune(cluster_list, bedpe, max_distance, min_cluster_size,print_ineligible,bedpe_out):
    new_cluster_list = []
    global cluster_counter
    for cluster in cluster_list:
        # cluster is beyond updatable window:
        if (bedpe is None   
            or cluster.chrom_a != bedpe.chrom_a
            or cluster.min_a - max_distance > bedpe.end_a
            or cluster.max_a + max_distance < bedpe.start_a):
            
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
        bedpe = Bedpe(line)
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

# --------------------------------------
# main function
def main():
    # parse the command line args
    args = get_args()
    # call primary function
    cluster_bedpe(args.input,
                  args.max_distance,
                  args.eval_param,
                  args.output,
                  args.is_sorted)

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
