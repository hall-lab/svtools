import argparse
import sys
from subprocess import Popen, PIPE, STDOUT
from svtools.vcf.file import Vcf
import svtools.utils as su

def run_cnvnator(cnvnator_path, root, window, coord_list):
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
    return cn_list

def sv_readdepth(vcf_file, sample, root, window, vcf_out, cnvnator_path, coord_list):
    cn = run_cnvnator(cnvnator_path, root, window, coord_list)
    write_copynumber(vcf_file, sample, vcf_out, cn)

def write_copynumber(vcf_file, sample, vcf_out, cn_list):
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
                        sys.stderr.write("Please input valid VCF, format field for " + sample + " not found in VCF")
                        sys.exit(1)
                  line = '\t'.join(map(str, line.rstrip().split('\t')[:9] + [sample]))
                  header.append(line)
                  continue
            else:
                in_header = False
                vcf.add_header(header)
                vcf.add_format('CN', 1, 'Float', 'Copy number of structural variant segment.')
                vcf_out.write(vcf.get_header() + '\n')
        v = line.rstrip().split('\t')
        # XXX Is this second check necessary? Wouldn't this be handled above? Missing header would hit this?
        if s_index == -1:
            sys.stderr.write("Input a valid sample name: " + sample + " not found in a provided VCF")
            sys.exit(1)
        v = v[:9] + [v[s_index]]
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
    return 'add copynumber information using cnvnator-multi'

def epilog():
    return '''As this program runs cnvnator-multi you must provide its location and must remember to have the ROOT package installed and properly configured. The input VCF file may be gzipped. If the input VCF file is omitted then the tool reads from stdin. Note that the coordinates file must end with a line containing the word exit.'''

def add_arguments_to_parser(parser):
    parser.add_argument('-c', '--coordinates', metavar='<FILE>', type=argparse.FileType('r'), required=True, default=None, help='file containing coordinate for which to retrieve copynumber (required)')
    parser.add_argument('-r', '--root', metavar='<FILE>', required=True, help='CNVnator .root histogram file (required)')
    parser.add_argument('-w', '--window', metavar='<INT>', required=True, help='CNVnator window size (required)')
    parser.add_argument('-s', '--sample', metavar='<STRING>', required=True, help='sample to annotate (required)')
    parser.add_argument('--cnvnator', metavar='<PATH>', required=True, help='path to cnvnator-multi binary (required)')
    parser.add_argument('-v', '--input-vcf', metavar='<VCF>', default=None, help='VCF input')
    parser.add_argument('-o', '--output-vcf', metavar='<PATH>', type=argparse.FileType('w'), default=sys.stdout, help='output VCF to write (default: stdout)')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description(), epilog=epilog())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.input_vcf) as stream:
        sv_readdepth(stream, args.sample, args.root, args.window, args.output_vcf, args.cnvnator, args.coordinates)

# initialize the script
if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
