import re
import sys

class Bedpe(object):
    def __init__(self, bed_list):
        self.c1 = bed_list[0]
        self.s1 = int(bed_list[1])
        self.e1 = int(bed_list[2])
        self.c2 = bed_list[3]
        self.s2 = int(bed_list[4])
        self.e2 = int(bed_list[5])
        self.name = bed_list[6]
        self.score = self.parse_score(bed_list[7])
        self.o1 = bed_list[8]
        self.o2 = bed_list[9]
        self.svtype = bed_list[10]
        self.filter = bed_list[11]
        self.malformedFlag = 0
        self.info1 = bed_list[12]
        self.info2 = bed_list[13]
        self.misc = bed_list[14:]
        self.check_malformed()

        # FIXME This is only really needed for varlookup. Something more general would be helpful
        self.cohort_vars = dict()
        
        try:
            self.svtype = self.retrieve_svtype()
        except ValueError:
            sys.stderr.write('No SVTYPE parseable for {0}'.format('\t'.join(bed_list)))
            sys.exit(1)
        self.af = self.retrieve_af()
        if self.svtype != bed_list[10]:
            sys.stderr.write("SVTYPE at Column 11({0})) and SVTYPE in INFO Column({1}) don't match at variant ID {3}\n".format(str(bed_list[10]), str(self.svtype), self.name))

    @staticmethod
    def parse_score(score):
        if score.isdigit():
            return float(score)
        else:
            return score

    def check_malformed(self):
        if self.info1 == 'MISSING':
            self.malformedFlag = 1
            self.info1 = self.info2
        if self.info2 == 'MISSING':
            self.malformedFlag = 2      

    def retrieve_svtype(self):
        try:
            svtype = re.split('=', ''.join(filter(lambda x: 'SVTYPE=' in x, self.info1.split(';'))))[1]
        except IndexError:
            raise ValueError
        return svtype

    def retrieve_af(self):
        try:
            af = re.split('=', ''.join(filter(lambda x: 'AF=' in x, self.info1.split(';'))))[1]
        except IndexError:
            af = None
        return af
    
    def __str__(self):
        '''
        A string representation of the line represented by this object
        '''
        return '\t'.join([
            self.c1,
            str(self.s1),
            str(self.e1),
            self.c2,
            str(self.s2),
            str(self.e2),
            self.name,
            str(self.score),
            self.o1,
            self.o2,
            self.svtype,
            self.filter,
            self.info1,
            self.info2] + 
            self.misc
            )

