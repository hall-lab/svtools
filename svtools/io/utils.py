import sys

class InputStream(object):
    '''This class handles opening either stdin or a gzipped or non-gzipped file'''
    def __init__(self, string):
        if string in (None, '-') and self.valid(string):
            self.handle = sys.stdin
        elif string.endswith('.gz'):
            import gzip
            self.handle = gzip.open(string, 'rb')
        else:
            self.handle = open(string, 'r')

    @staticmethod
    def valid(string):
        if string is None and sys.stdin.isatty():
            raise IOError('no input specified but terminal is interactive')
        return True
    
    def __enter__(self):
        return self.handle

    def __exit__(self, *kwargs):
        self.handle.close()

    def __iter__(self):
        return self.handle.__iter__()

