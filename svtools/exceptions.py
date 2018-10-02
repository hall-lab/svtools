class SvtoolsException(Exception):
    pass

class MissingProbabilitiesException(SvtoolsException):
    def __init__(self, string):
        msg  = (
                "Please ensure you've run lumpy with"
                "the -P option to emit breakpoint probabilities."
                )
        super(MissingProbabilitiesException, self).__init__(
                ' '.join((string, msg))
                )
