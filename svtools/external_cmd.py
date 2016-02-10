from abc import ABCMeta
import sys
import os
import argparse
import subprocess
import signal

class ShellCmd(object):
    __metaclass__ = ABCMeta

    def __init__(self, name, relative_path=None):
        self.name = name
        if relative_path is not None:
            self.relative_path = relative_path
        else:
            self.relative_path = name

    def path_to_shell_script(self):
        # FIXME This may not be the best way to find the location of these scripts. See pkg_resources as a possible alternative
        path = os.path.dirname(os.path.abspath(__file__))

        path_to_script = os.path.join(path, self.relative_path)
        if os.path.isfile(path_to_script):
            return path_to_script
        else:
            sys.stderr.write('Unable to locate {0} script at {1}\n'.format(self.name, path_to_script))
            sys.exit(1)
    
    def run_cmd_with_options(self, options):
        cmd = [ self.path_to_shell_script() ]
        cmd.extend(options)
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
                sys.stderr.write('{0} exited with code {1}\n'.format(self.name, code))
            sys.exit(code)
