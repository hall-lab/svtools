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

    @staticmethod
    def parse_info_tag(info_string, tag):
        '''
        Accessory method to parse out the value of a tag in an info string.
        Make sure to include the equals sign if you are looking for a
        non-boolean tag
        '''

        tag_start = info_string.find(tag)
        if tag_start == -1:
            # If you were looking for a flag then this is the right value.
            # Otherwise your tag doesn't exist. Client code must know how to
            # interpret.
            return False

        tag_end = info_string.find(';', tag_start)
        value_start = tag_start + len(tag)
        if (value_start >= len(info_string)) or (tag_end != -1 and value_start >= tag_end):
            return True
        if tag_end == -1:
            tag_end = None # We didn't find our end index
        return info_string[value_start:tag_end]

    @staticmethod
    def update_info_tag(info_string, tag, new_value):
        '''
        Accessory method to update a tag's value. Like parse_info_tag, make sure to include the equals sign.
        '''

        tag_start = info_string.find(tag)
        if tag_start == -1:
            raise ValueError("Tag {0} doesn't exist".format(tag))

        tag_end = info_string.find(';', tag_start)
        value_start = tag_start + len(tag)
        if (value_start >= tag_end and tag_end != -1) or value_start >= len(info_string):
            raise ValueError("Tag {0} doesn't have a value".format(tag))
        if tag_end == -1:
            tag_end = None # We didn't find our end index
        new_info_string = info_string[:value_start] + new_value
        if tag_end:
            new_info_string += info_string[tag_end:]
        return new_info_string

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
            raise ValueError('SVTYPE field not present in INFO field')
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

    @staticmethod
    def sname_value(info_string):
        '''
        Retrieves the SNAME value from an info_string. Static method so we can
        easily do it for either info1 or info2 on demand
        '''
        value = Bedpe.parse_info_tag(info_string, 'SNAME=')
        if value in (False, True):
            return None
        else:
            return value

    @staticmethod
    def _combine_sname_values(first, second):
        '''
        Combine the sname values from two comma-separated strings
        '''
        combined = None
        if first is not None and second is not None:
            sname_set = set(first.split(',') + second.split(','))
            combined = ','.join(sname_set)
        else:
            combined = first or second # set to whichever is non-None
        return combined

    @staticmethod
    def _update_sname_field(original_info1, original_info2):
        '''
        Update the sname field in the original info string by adding
        values from the another info string
        '''
        new_sname = Bedpe._combine_sname_values(
                Bedpe.sname_value(original_info1),
                Bedpe.sname_value(original_info2))
        if new_sname:
            return Bedpe.update_info_tag(original_info1, 'SNAME=', new_sname)
        else:
            return original_info1

    def combine_snames(self, other):
        '''
        This method adds the sname values from the info fields of another bedpe
        entry into the sname tag of itself
        '''
        self.info1 = self._update_sname_field(self.info1, other.info1)
        self.info2 = self._update_sname_field(self.info2, other.info2)
