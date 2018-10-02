import l_bp
from exceptions import MissingProbabilitiesException

class BreakpointInterval(object):
    '''
    Class for storing the range and probability distribution
    of a breakpoint
    '''
    # Constant value for slop padding
    SLOP_PROB = 1e-100

    def __init__(self, chrom, start, end, p):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.p = p

    def pad_slop(self, percent_slop, fixed_slop):
        '''
        Add slop to the interval
        '''
        slop = int(max(percent_slop * (self.end - self.start + 1), fixed_slop))
        self.start -= slop
        self.end += slop
        self.p = [BreakpointInterval.SLOP_PROB] * slop + self.p + [BreakpointInterval.SLOP_PROB] * slop
        self._trim()
        self._normalize()

    def _trim(self):
        '''
        Trim any part of range past the beginning of the chromosome
        '''
        if self.start < 0:
            self.p = self.p[-self.start:]
            self.start = 0

    def _normalize(self):
        '''
        Normalize interval's probability to sum to 1
        '''
        sum_p = sum(self.p)
        self.p = [float(x)/sum_p for x in self.p]

    def common_range(self, other):
        return max(self.start, other.start), min(self.end, other.end)

    def overlap_prob(self, other, c_start, c_len):
        start_off = c_start - self.start
        other_start_off = c_start - other.start
        ovl = 0
        for i in range(c_len):
            ovl += min(self.p[i + start_off], other.p[i + other_start_off])
        return ovl


class Breakpoint(object):
    '''
    Class for storing information about Breakpoints for merging
    '''

    def __init__(self, line, percent_slop=0, fixed_slop=0):
        '''
        Initialize with slop for probabilities
        '''
        self.l = line

        (self.sv_type,
        chr_l,
        chr_r,
        self.strands,
        start_l,
        end_l,
        start_r,
        end_r,
        m) = l_bp.split_v(line)

        try:
            self.left = BreakpointInterval(chr_l, start_l, end_l, self.floats_from_tag(m, 'PRPOS'))
            self.right = BreakpointInterval(chr_r, start_r, end_r, self.floats_from_tag(m, 'PREND'))
        except RuntimeError as e:
            raise MissingProbabilitiesException(str(e))

        if ((percent_slop > 0) or (fixed_slop > 0)):
            self.left.pad_slop(percent_slop, fixed_slop)
            self.right.pad_slop(percent_slop, fixed_slop)

    def __str__(self):
        '''
        Convert back to a string
        '''
        return '\t'.join([str(x) for x in [self.left.chrom,
                                           self.left.start,
                                           self.left.end,
                                           self.right.chrom,
                                           self.right.start,
                                           self.right.end,
                                           self.sv_type,
                                           self.strands,
                                           self.left.p,
                                           self.right.p]])
    def ovl(self, b):
        '''
        Calculate overlapping cumulative probability value as weight?
        0 if not overlapping.
        '''
        if ((self.left.chrom != b.left.chrom) or
            (self.right.chrom != b.right.chrom) or
            (self.sv_type != b.sv_type)):
                return 0
        #get common intervals
        c_start_l, c_end_l = self.left.common_range(b.left)
        c_start_r, c_end_r = self.right.common_range(b.right)

        c_l_len = c_end_l - c_start_l + 1
        c_r_len = c_end_r - c_start_r + 1

        if (c_l_len < 1) or (c_r_len < 1):
            return 0

        ovl_l = self.left.overlap_prob(b.left, c_start_l, c_l_len)
        ovl_r = self.right.overlap_prob(b.right, c_start_r, c_r_len)

        return ovl_l * ovl_r

    @staticmethod
    def floats_from_tag(info_dict, tag):
        if tag in info_dict:
            return [float(x) for x in info_dict[tag].split(',')]
        else:
            raise RuntimeError('Required tag {0} not found.'.format(tag))

