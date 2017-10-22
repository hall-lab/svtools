import l_bp

class Breakpoint:
    '''
    Class for storing information about Breakpoints for merging
    '''
    # Constant value for slop padding
    SLOP_PROB = 1e-100

    def __init__(self, line, percent_slop=0, fixed_slop=0):
        '''
        Initialize with slop for probabilities
        '''
        self.l = line

        (self.sv_type,
        self.chr_l,
        self.chr_r,
        self.strands,
        self.start_l,
        self.end_l,
        self.start_r,
        self.end_r,
        m) = l_bp.split_v(line)

        self.p_l = self.floats_from_tag(m, 'PRPOS')
        self.p_r = self.floats_from_tag(m, 'PREND')

        if ((percent_slop > 0) or (fixed_slop > 0)):

            self.start_l, self.end_l, new_p_l = self.pad_slop(self.start_l, self.end_l, self.p_l, percent_slop, fixed_slop)
            self.start_r, self.end_r, new_p_r = self.pad_slop(self.start_r, self.end_r, self.p_r, percent_slop, fixed_slop)

            # chew off overhang if self.start_l or self.start_r less than 0
            self.start_l, new_p_l = self.trim_slop(self.start_l, new_p_l)
            self.start_r, new_p_r = self.trim_slop(self.start_r, new_p_r)

            # normalize so each probability curve sums to 1.
            self.p_l = self.normalize(new_p_l)
            self.p_r = self.normalize(new_p_r)

    def __str__(self):
        '''
        Convert back to a string
        '''
        return '\t'.join([str(x) for x in [self.chr_l,
                                           self.start_l,
                                           self.end_l,
                                           self.chr_r,
                                           self.start_r,
                                           self.end_r,
                                           self.sv_type,
                                           self.strands,
                                           self.p_l,
                                           self.p_r]])
    def ovl(self, b):
        '''
        Calculate overlapping cumulative probability value as weight?
        0 if not overlapping.
        '''
        if ((self.chr_l != b.chr_l) or
            (self.chr_r != b.chr_r) or
            (self.sv_type != b.sv_type)):
                return 0
        #get common intervals
        c_start_l, c_end_l = self.common_interval(self.start_l, b.start_l, self.end_l, b.end_l)
        c_start_r, c_end_r = self.common_interval(self.start_r, b.start_r, self.end_r, b.end_r)

        c_l_len = c_end_l - c_start_l + 1
        c_r_len = c_end_r - c_start_r + 1

        if (c_l_len < 1) or (c_r_len < 1):
            return 0

        ovl_l = self.overlap_prob(c_l_len, c_start_l, self.start_l, b.start_l, self.p_l, b.p_l)
        ovl_r = self.overlap_prob(c_r_len, c_start_r, self.start_r, b.start_r, self.p_r, b.p_r)

        return ovl_l * ovl_r

    @staticmethod
    def common_interval(start1, start2, end1, end2):
        return max(start1, start2), min(end1, end2)

    @staticmethod
    def overlap_prob(c_len, c_start, start, b_start, p, b_p):
        start_off = c_start - start
        b_start_off = c_start - b_start
        ovl = 0
        for i in range(c_len):
            ovl += min(p[i + start_off], b_p[i + b_start_off])
        return ovl

    @staticmethod
    def floats_from_tag(info_dict, tag):
        if tag in info_dict:
            return [float(x) for x in info_dict[tag].split(',')]
        else:
            raise RuntimeError('Required tag {0} not found.'.format(tag))

    @staticmethod
    def pad_slop(start, end, p, percent_slop, fixed_slop):
        slop = int(max(percent_slop * (end - start + 1), fixed_slop))
        new_start = start - slop
        new_end = end + slop
        new_p = [Breakpoint.SLOP_PROB] * slop + p + [Breakpoint.SLOP_PROB] * slop
        return (new_start, new_end, new_p)

    @staticmethod
    def trim_slop(start, p):
        if start < 0:
            new_p = p[-start:]
            return (0, new_p)
        else:
            return (start, p)

    @staticmethod
    def normalize(p):
        sum_p = sum(p)
        return [float(x)/sum_p for x in p]
