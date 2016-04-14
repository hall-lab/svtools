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
        self.orig_name1 = bed_list[12]
        self.orig_ref1 = bed_list[13]
        self.orig_alt1 = bed_list[14]
        self.orig_name2 = bed_list[15]
        self.orig_ref2 = bed_list[16]
        self.orig_alt2 = bed_list[17]
        self.malformedFlag = 0
        self.info1 = bed_list[18]
        self.info2 = bed_list[19]
        self.misc = bed_list[20:]
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
    @property
    def info(self):
        '''
        Return the appropriate info field if only one is required. Info from the primary variant is preferred if available.
        '''
        if self.info1 == 'MISSING':
            return self.info2
        else:
            return self.info1

    def set_info(self, field, value):
        '''
        Add the info field to the BEDPE line info fields. As BEDPE lines don't know about their headers this is not a safe operation.
        Doesn't add to info field if it is the null character. Probably this is wrong.
        '''
        new_tag = ';' + str(field);
        if value is not None:
            new_tag += '=' + str(value)
        if self.malformedFlag != 1:
            self.info1 = self.info1 + new_tag
        if self.malformedFlag != 2 and self.info2 != '.':
            self.info2 = self.info2 + new_tag

    def check_malformed(self):
        if self.info1 == 'MISSING':
            self.malformedFlag = 1
        if self.info2 == 'MISSING':
            self.malformedFlag = 2      

    def retrieve_svtype(self):
        try:
            svtype = re.split('=', ''.join(filter(lambda x: 'SVTYPE=' in x, self.info.split(';'))))[1]
        except IndexError:
            raise ValueError
        return svtype

    def retrieve_af(self):
        try:
            af = re.split('=', ''.join(filter(lambda x: 'AF=' in x, self.info.split(';'))))[1]
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
            self.orig_name1,
            self.orig_ref1,
            self.orig_alt1,
            self.orig_name2,
            self.orig_ref2,
            self.orig_alt2,
            self.info1,
            self.info2] + 
            self.misc
            )

