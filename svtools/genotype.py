import sys
import os
import argparse
import subprocess
import signal
from collections import namedtuple

def path_to_shell_script():
    # FIXME This may not be the best way to find the location of these scripts. See pkg_resources as a possible alternative
    path = os.path.dirname(os.path.abspath(__file__))
    path_to_script = os.path.join(path, 'bin', 'svtyper', 'svtyper')
    if os.path.isfile(path_to_script):
        return path_to_script
    else:
        sys.stderr.write("Unable to locate svtyper script")
        sys.exit(1)

def svtyper_option_lut():
    opts = { 
            'bam' : '-B',
            'split_bam' : '-S',
            'input_vcf' : '-i',
            'output_vcf' : '-o',
            'splflank' : '-f',
            'discflank' : '-F',
            'split_weight' : '--split_weight',
            'disc_weight' : '--disc_weight',
            'num_samp' : '-n',
            'legacy' : '-M',
            'debug' : '--debug',
            }
    return opts

def description():
    return 'compute genotype of structural variants based on breakpoint depth'

def add_arguments_to_parser(parser):
    parser.add_argument('-B', '--bam', type=str, required=True, help='BAM file(s), comma-separated if genotyping multiple BAMs')
    parser.add_argument('-S', '--split_bam', type=str, required=False, help='split-read bam file for sample, comma-separated if genotyping multiple BAMs')
    parser.add_argument('-i', '--input_vcf', help='VCF input (default: stdin)')
    parser.add_argument('-o', '--output_vcf', help='output VCF to write (default: stdout)')
    parser.add_argument('-f', '--splflank', type=int, required=False, default=20, help='min number of split read query bases flanking breakpoint on either side [20]')
    parser.add_argument('-F', '--discflank', type=int, required=False, default=20, help='min number of discordant read query bases flanking breakpoint on either side. (should not exceed read length) [20]')
    parser.add_argument('--split_weight', type=float, required=False, default=1, help='weight for split reads [1]')
    parser.add_argument('--disc_weight', type=float, required=False, default=1, help='weight for discordant paired-end reads [1]')
    parser.add_argument('-n', dest='num_samp', type=int, required=False, default=1000000, help='number of pairs to sample from BAM file for building insert size distribution [1000000]')
    parser.add_argument('-M', action='store_true', dest='legacy', required=False, help='split reads are flagged as secondary, not supplementary. For compatibility with legacy BWA-MEM "-M" flag')
    parser.add_argument('--debug', action='store_true', help='debugging verbosity')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    cmd = [ path_to_shell_script() ]
    argdict = vars(args)
    optlut = svtyper_option_lut()
    for variable, value in argdict.iteritems():
        if variable != 'entry_point' and value not in (False, None):
            cmd.extend([optlut[variable], str(value)])
    # Here we re-set Python's treatment of SIGPIPE to the default
    # as described here: http://www.chiark.greenend.org.uk/~cjwatson/blog/python-sigpipe.html
    sys.stderr.write('Running svtyper with options: {0}\n'.format(' '.join(cmd)))
    p = subprocess.Popen(cmd, preexec_fn=lambda:
            signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    code = p.wait()
    # The check for code != 141 is here because 
    # 141 indicates a SIGPIPE signal returned in the underlying bash pipelines
    # We want to be silent there
    # FIXME 141 is bash specific and while the underlying scripts are bash
    # It is not clear that the shell this script is run in should be exiting with 141 
    # or if that is even necessary
    if code:
        if code is not 141:
            sys.stderr.write('svtyper exited with code {0}\n'.format(code))
        sys.exit(code)

if __name__ == "__main__":
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
