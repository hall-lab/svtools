import sys

class InputStream(object):
    '''This class handles opening either stdin or a gzipped or non-gzipped file'''

    def __init__(self, string):
        '''Create a new wrapper around a stream'''
        if string in (None, '-', 'stdin') and self.valid(string):
            self.handle = sys.stdin
        elif string.endswith('.gz'):
            import gzip
            self.handle = gzip.open(string, 'rb')
        else:
            self.handle = open(string, 'r')

    @staticmethod
    def valid(string):
        '''Check if we allow the opening of stdin'''
        if string is None and sys.stdin.isatty():
            raise IOError('no input specified but terminal is interactive')
        return True
    
    def __enter__(self):
        '''Support use of with by passing back the originating handle'''
        return self.handle

    def __exit__(self, *kwargs):
        '''Support use of with by closing on exit of the context'''
        self.handle.close()

    def __iter__(self):
        '''Support use in loops like a normal file object'''
        return self.handle.__iter__()

    def close(self):
        '''Close the underlying handle'''
        return self.handle.close()
