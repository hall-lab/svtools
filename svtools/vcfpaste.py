import argparse, sys
import gzip

MAX_SPLIT = 9

class Vcfpaste(object):
    def __init__(self, vcf_list, master=None, sum_quals=None):
        self.vcf_list = vcf_list
        self.master = master
        self.sum_quals = sum_quals

    def execute(self, output_handle=sys.stdout):
        try:
            self.read_filenames()
            self.open_files()
            self.write_header(output_handle)
            self.write_variants(output_handle)
        finally:
            self.close_files()

    def read_filenames(self):
        self.vcf_file_names = []
        with open(self.vcf_list, 'r') as vcf_list_file:
            for line in vcf_list_file:
                path = line.rstrip()
                self.vcf_file_names.append(path)
        if self.master == None:
            self.master = self.vcf_file_names[0]
        self.vcf_file_names.insert(0, self.master)

    def open_files(self):
        self.vcf_files = []
        # parse the vcf files to paste
        for path in self.vcf_file_names:
            if path.endswith('.gz'):
                self.vcf_files.append(gzip.open(path, 'rb'))
            else:
                self.vcf_files.append(open(path, 'r'))
    
    def write_header(self, output_handle=sys.stdout):
        master = self.vcf_files[0]
        while 1:
            master_line = master.readline()
            if not master_line:
                break
            if master_line[:2] != '##':
                break
            output_handle.write(master_line)
        out_v = master_line.rstrip().split('\t', MAX_SPLIT)[:(MAX_SPLIT - 1)] + ["FORMAT"]

        for vcf in self.vcf_files[1:]:
            while 1:
                l = vcf.readline()
                if not l:
                    break
                if l[:2] == '##':
                    continue
                if l[0] == '#':
                    out_v = out_v + l.rstrip().split('\t', MAX_SPLIT)[MAX_SPLIT:]
                    break
        output_handle.write('\t'.join(map(str, out_v)) + '\n')

    def write_variants(self, output_handle=sys.stdout):
        while 1:
            master_line = self.vcf_files[0].readline()
            if not master_line:
                break
            master_v = master_line.rstrip().split('\t', MAX_SPLIT)
            if len(master_v) < 8:
                sys.stderr.write('\nERROR: Master file {0} had less than 8 columns.\n'.format(self.vcf_files[0].name))
                exit(1)
            out_v = master_v[:8] # output array of fields
            qual = 0
            if out_v[5] != '.':
                qual = float(out_v[5])
            format = None # column 9, VCF format field.

            for vcf in self.vcf_files[1:]:
                line = vcf.readline()
                if not line:
                    sys.stderr.write('\nERROR: VCF files differ in length\n')
                    exit(1)
                line_v = line.rstrip().split('\t', MAX_SPLIT)
                if len(line_v) < 10:
                    sys.stderr.write('\nERROR: {0} had less than 10 columns. Only the master may be an 8 column VCF.\n'.format(vcf.name))
                    exit(1)

                # set FORMAT field as format in first VCF.
                # cannot extract this from master, since it may have
                # been altered in the processing of the VCFs.
                if format is None:
                    format = line_v[8]
                    out_v.append(format)

                if line_v[5] != '.':
                    qual += float(line_v[5])
                out_v = out_v + line_v[9:]
            if self.sum_quals:
                out_v[5] = qual
            output_handle.write( '\t'.join(map(str, out_v)) + '\n')

    def close_files(self):
        for f in self.vcf_files:
            f.close()

def description():
    return 'paste VCFs from multiple samples'

def epilog():
    return '''VCF files may be gzipped. If the -m argument is omitted then the first file in the list of files in --vcf-list is treated as the master.'''

def add_arguments_to_parser(parser):
    parser.add_argument('-f', '--vcf-list', metavar='<FILE>', required=True, help='file containing a line-delimited list of VCF files to paste (required)')
    parser.add_argument('-m', '--master', metavar='<VCF>', default=None, help='VCF file to set first 8 columns of variant info (otherwise first file in --vcf-list)')
    parser.add_argument('-q', '--sum-quals', required=False, action='store_true', help='sum QUAL scores of input VCFs as output QUAL score')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description(), epilog=epilog())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    paster = Vcfpaste(args.vcf_list, master=args.master, sum_quals=args.sum_quals)
    paster.execute()

# initialize the script
if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
