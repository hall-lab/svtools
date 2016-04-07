import argparse, sys
from svtools.bedpe import Bedpe
from svtools.cluster import Cluster
import svtools.utils as su

class Pruner(object):
    def __init__(self, max_distance, eval_param):
        self.max_distance = max_distance
        self.eval_param = eval_param
        self.min_cluster_size = 1
        self.bedpe_lines = 0
        self.skipped_lines = 0
        self.emitted_lines = 0
        self.cluster_list = []

    def stats_report_string(self):
        return "{0} variants read\n{1} skipped\n{2} pruned\n{3} emitted\n".format(self.bedpe_lines,
            self.skipped_lines,
            self.bedpe_lines - self.skipped_lines - self.emitted_lines,
            self.emitted_lines)

    def cluster_bedpe(self, in_file, bedpe_out, is_sorted):
        # Locally alias instance variables
        max_distance = self.max_distance
        eval_param = self.eval_param

        in_header = True
        for line in in_file:
            if line.startswith('#') and in_header:
                bedpe_out.write(line)
                continue
            in_header = False
            self.bedpe_lines += 1
            bedpe = Bedpe(line.rstrip().split('\t'))
            if bedpe.af is None:
                sys.stderr.write('No allele frequency for variant found. This tool requires allele frequency information to function. Please add with svtools afreq and rerun\n')
                sys.exit(1)
            if bedpe.af == '.':
                self.skipped_lines += 1
                continue
            matched_clusters = []
            for cluster in self.cluster_list:
                if cluster.can_add(bedpe, max_distance):
                    cluster.add(bedpe, eval_param)
                    matched_clusters.append(cluster)
            if not matched_clusters:
                new_cluster = Cluster()
                new_cluster.add(bedpe, eval_param)
                self.cluster_list.append(new_cluster)
            else:
                if len(matched_clusters) > 1:
                    i = 0
                    matched_cluster_pruned = False
                    while i < (len(matched_clusters) - 1):
                        j = i + 1
                        while j < len(matched_clusters):
                            if matched_clusters[i].can_add(matched_clusters[j].elements[0], max_distance):
                                matched_clusters[i].add(matched_clusters[j].elements[0], eval_param)
                                matched_cluster_pruned = True
                                del matched_clusters[j]
                            j += 1
                        i += 1        
                    if matched_cluster_pruned:
                        self.cluster_list = [cluster for cluster in self.cluster_list if cluster not in matched_clusters]
            #prune and print eligible clusters
            if self.bedpe_lines % 1000 == 0 and is_sorted:
                self.cluster_list = self.prune(bedpe,
                                     False,
                                     bedpe_out)
    
        self.cluster_list = self.prune(None,
                             True,
                             bedpe_out)

        sys.stderr.write(self.stats_report_string())
        return

    def prune(self, bedpe, print_ineligible, bedpe_out):
        '''
        Prints and removes clusters from the cluster_list that are beyond the distance window
        '''
        min_cluster_size = self.min_cluster_size
        max_distance = self.max_distance
        new_cluster_list = []
        for cluster in self.cluster_list:
            # cluster is beyond updatable window:
            if (bedpe is None or 
                    cluster.chrom_a != bedpe.c1 or 
                    cluster.min_a - max_distance > bedpe.e1 or 
                    cluster.max_a + max_distance < bedpe.s1):
                
                # print the cluster if eligible
                if (cluster.size >= min_cluster_size or print_ineligible):
                    self.emitted_lines += 1
                    bedpe_out.write(cluster.get_cluster_string() + '\n')
                else:
                    new_cluster_list.append(cluster)
            # cluster is still within updatable window,
            # leave it in the cluster list
            else:
                new_cluster_list.append(cluster)
    
        return new_cluster_list

def description():
    return 'cluster and prune a BEDPE file by position based on allele frequency'

def epilog():
    return 'The input BEDPE file can be gzipped if it is specified explicitly.'

def add_arguments_to_parser(parser):
    parser.add_argument('-d', '--distance', metavar='<INT>', type=int, dest='max_distance', default=50, help='max separation distance (bp) of adjacent loci in cluster [50]')
    parser.add_argument('-e', '--eval-param', metavar='<STRING>', help='evaluating parameter for choosing best bedpe in a cluster(e.g. af=AlleleFrequency default:af)')
    parser.add_argument('-s', '--is-sorted', action='store_true', help='specify if an input file is sorted. Sort with svtools bedpesort. (default=False)')
    parser.add_argument('input', nargs='?', metavar='<BEDPE>', default=None, help='BEDPE file to read. If \'-\' or absent then defaults to stdin.')
    parser.add_argument('-o', '--output', metavar='<BEDPE>', type=argparse.FileType('w'), default=sys.stdout, help='output bedpe to write (default: stdout)')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description(), epilog=epilog())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.input) as stream:
        pruner = Pruner(args.max_distance, args.eval_param)
        pruner.cluster_bedpe(stream, args.output, args.is_sorted)

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
