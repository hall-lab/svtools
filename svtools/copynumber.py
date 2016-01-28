import argparse
import sys
from subprocess import Popen, PIPE, STDOUT
from svtools.vcf.file import Vcf


def sv_readdepth(vcf_file, sample, root, window, vcf_out, cnvnator_path, coord_list):
    # Read and run cnvnator
    p1 = Popen(['cat', coord_list.name], stdout=PIPE)
    cmd = map(str, [cnvnator_path, '-root', root, '-genotype', window])
    p2 = Popen(cmd, stdin=p1.stdout, stdout=PIPE)
    # NOTE the below awk command selects only those lines which were genotyped
    # The line beginning 'Assuming' looks to be the program reporting its assuming a male
    # Each line looks like
    # Genotype chr1:1-13000 5173T.root 0.377524 0.368896
    # The first number is the copy number for the requested window
    # The second seems to be the copy number for a fixed window size of 1000
    # See http://wiki.biouml.org/index.php/CNVnator_genotype_output_(file_format)
    p3 = Popen(['awk', '{ if($1!="Assuming"){print $4} }'], stdin=p2.stdout, stdout=PIPE)
    cn_list = map(float, p3.communicate()[0].split('\n')[:-1])

    #go through the VCF and add the read depth annotations
    in_header = True
    header = []
    vcf = Vcf()
    i = 0
    s_index = -1
    for line in vcf_file:
        if in_header:
            if line[0] == '#' and line[1] == '#':
                header.append(line)
                continue
            if line[0] == '#' and line[1] != '#':
                  try:
                        s_index = line.rstrip().split('\t').index(sample)
                  except ValueError:
                        stderr.write("Please input valid VCF, format field for " + sample + " not found in VCF")
                        sys.exit(1)
                  line = '\t'.join(map(str,[line.rstrip().split('\t')[0:8],sample]))
                  header.append(line)
                  continue
            else:
                in_header = False
                vcf.add_header(header)
                vcf.add_format('CN', 1, 'Float', 'Copy number of structural variant segment.')
                vcf_out.write('\t'.join([vcf.get_header(include_samples=False), 'FORMAT', sample]) + '\n')
        v = line.rstrip().split('\t')
        # XXX Is this second check necessary? Wouldn't this be handled above? Missing header would hit this?
        if s_index == -1:
            stderr.write("Input a valid sample name: " + sample + " not found in a provided VCF")
            sys.exit(1)
        v = v[0:9] + [v[s_index]]
        if not any("SVTYPE=BND" in s for s in v):
            if "CN" not in v[8]:
                v[8] = v[8] + ":CN"
                v[9] = v[9] + ":" + str(cn_list[i])
            else:
                cn_index = v[8].rstrip().split(":").index("CN")
                gts = v[9].rstrip().split(":")
                gts[cn_index] = str(cn_list[i])
                v[9] = ":".join(gts)
            i += 1
        # write the VCF
        vcf_out.write('\t'.join(v) + '\n')
    vcf_out.close()
    return

def description():
    return 'Compute genotype of structural variants based on breakpoint depth'

def add_arguments_to_parser(parser):
    parser.add_argument('-v', '--input-vcf', type=argparse.FileType('r'), default=None, help='VCF input')
    parser.add_argument('-c', '--coordinates', type=argparse.FileType('r'), required=True, default=None, help='BED input')
    parser.add_argument('-r', '--root', required=True, help='CNVnator .root histogram file (required)')
    parser.add_argument('-w', '--window', required=True, help='CNVnator window size (required)')
    parser.add_argument('-s', '--sample', required=True, help='sample to annotate')
    parser.add_argument('--cnvnator', required=True, help='path to cnvnator-multi binary')
    parser.add_argument('-o', '--output-vcf', type=argparse.FileType('w'), default=sys.stdout, help='output VCF to write (default: stdout)')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    if args.input_vcf == None:
        if sys.stdin.isatty():
            parser.print_help()
            sys.exit(1)
        else:
            args.input_vcf = sys.stdin
    sv_readdepth(args.input_vcf, args.sample, args.root, args.window, args.output_vcf, args.cnvnator, args.coordinates)
    if args.input_vcf != sys.stdin:
        args.input_vcf.close()

# initialize the script
if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))