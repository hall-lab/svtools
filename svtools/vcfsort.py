import sys
import os
import argparse
import subprocess
import signal

def path_to_shell_script():
    # FIXME This may not be the best way to find the location of these scripts. See pkg_resources as a possible alternative
    # or possible try this
    #path = os.path.dirname(sys.modules['svtools'].__file__)
    #path_to_script = os.path.join(path, 'bin', 'vcfsort')
    path = os.path.dirname(os.path.abspath(__file__))
    path_to_script = os.path.join(path, 'bin', 'vcfsort')
    if os.path.isfile(path_to_script):
        return path_to_script
    else:
        sys.write.stderr("Unable to locate vcf script")
        sys.exit(1)

def description():
    return 'sort a VCF file'

def add_arguments_to_parser(parser):
    parser.add_argument('input', metavar='<VCF file>', nargs='?', help='VCF file to sort')
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
    # Here we re-set Python's treatment of SIGPIPE to the default
    # as described here: http://www.chiark.greenend.org.uk/~cjwatson/blog/python-sigpipe.html
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
            sys.stderr.write('vcfsort bash script exited with code {0}\n'.format(code))
        sys.exit(code)

if __name__ == "__main__":
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
