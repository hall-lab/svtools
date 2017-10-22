import l_bp

class Breakpoint:
    '''
    Class for storing information about Breakpoints for merging
    '''
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

        slop_prob = 1e-100 # FIXME This is a constant. Pull out to make more obvious
        if ((percent_slop > 0) or (fixed_slop > 0)):

            l_slop = int(max(percent_slop * (self.end_l - self.start_l + 1), fixed_slop))
            r_slop = int(max(percent_slop * (self.end_r - self.start_r + 1), fixed_slop))

            # pad each interval with slop_prob on each side. TODO This should be a method
            self.start_l = self.start_l - l_slop
            self.end_l = self.end_l + l_slop
            new_p_l = [slop_prob] * l_slop + self.p_l + [slop_prob] * l_slop

            self.start_r = self.start_r - r_slop
            self.end_r = self.end_r + r_slop
            new_p_r = [slop_prob] * r_slop + self.p_r + [slop_prob] * r_slop

            # chew off overhang if self.start_l or self.start_r less than 0 TODO This should also be a method
            if self.start_l < 0:
                new_p_l = new_p_l[-self.start_l:]
                self.start_l = 0
            if self.start_r < 0:
                new_p_r = new_p_r[-self.start_r:]
                self.start_r = 0

            # normalize so each probability curve sums to 1. TODO Should be a method
            sum_p_l = sum(new_p_l)
            self.p_l = [float(x)/sum_p_l for x in new_p_l]
            sum_p_r = sum(new_p_r)
            self.p_r = [float(x)/sum_p_r for x in new_p_r]

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
        #get left common interval
        c_start_l, c_end_l = max(self.start_l, b.start_l), min(self.end_l, b.end_l)
        #get right common interval
        c_start_r, c_end_r = max(self.start_r, b.start_r), min(self.end_r, b.end_r)

        c_l_len = c_end_l - c_start_l + 1
        c_r_len = c_end_r - c_start_r + 1

        if (c_l_len < 1) or (c_r_len < 1):
            return 0

        # TODO This should probably be a method as well
        self_start_off_l = c_start_l - self.start_l
        b_start_off_l = c_start_l - b.start_l

        self_start_off_r = c_start_r - self.start_r
        b_start_off_r = c_start_r - b.start_r

        ovl_l = 0
        for i in range(c_l_len):
            ovl_l += min(self.p_l[i + self_start_off_l], b.p_l[i + b_start_off_l])

        ovl_r = 0
        for i in range(c_r_len):
            ovl_r += min(self.p_r[i + self_start_off_r], b.p_r[i + b_start_off_r])

        return ovl_l * ovl_r

    @staticmethod
    def floats_from_tag(info_dict, tag):
        if tag in info_dict:
            return [float(x) for x in info_dict[tag].split(',')]
        else:
            raise RuntimeError('Required tag {0} not found.'.format(tag))
