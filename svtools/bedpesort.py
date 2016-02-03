import sys
import os
import argparse
import subprocess
import signal

def path_to_shell_script():
    # FIXME This may not be the best way to find the location of these scripts. See pkg_resources as a possible alternative
    path = os.path.dirname(os.path.abspath(__file__))
    path_to_script = os.path.join(path, 'bin', 'bedpesort')
    if os.path.isfile(path_to_script):
        return path_to_script
    else:
        sys.write.stderr("Unable to locate bedpesort script")
        sys.exit(1)

def description():
    return 'sort a VCF file'

def add_arguments_to_parser(parser):
    parser.add_argument('input', metavar='<BEDPE file>', nargs='?', help='BEDPE file to sort')
    parser.add_argument('output', metavar='<output file>', nargs='?', help='output file to write to')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    cmd = [ path_to_shell_script() ]
    if args.input:
        cmd.append(args.input)
    if args.output:
        cmd.append(args.output)
    p = subprocess.Popen(cmd, preexec_fn=lambda:
            signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    code = p.wait()
    if code:
        sys.stderr.write('bedpesort bash script exited with code {0}\n'.format(code))
        sys.exit(code)

if __name__ == "__main__":
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
