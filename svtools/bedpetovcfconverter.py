from svtools.vcf.variant import Variant
import re

class BedpeToVcfConverter(object):
    '''
    This is a class to take Bedpe objects and convert them to VCF variant object(s)
    '''

    def __init__(self, vcf):
        '''
        Initialize a new converter with a reference to Vcf object
        '''
        self.vcf_header = vcf

    @staticmethod
    def adjust_by_tag(bedpe, info_tag, strand, coordinate):
        '''
        Undo adjustments to BEDPE coordinates based on strand and INFO tag values
        '''
        if strand == '-' and bedpe.svtype == 'BND':
            coordinate += 1
        if info_tag in bedpe.info:
            value_array = re.split('=|,', ''.join(filter(lambda x: info_tag in x, bedpe.info.split(';'))))
            coordinate -= int(value_array[1])
        return coordinate

    @staticmethod
    def determine_sep(o2):
        '''
        Given the orientation of the other breakend, determine separator for VCF alt
        '''
        if o2 == '+':
            return ']'
        else:
            return '['

    @staticmethod
    def determine_flanks(o1):
        '''
        Given the orientation of the breakend, return proper flanking sequence strings
        '''
        if o1 == '+':
            return ('N', '')
        else:
            return ('', 'N')

    def bnd_alt_string(self, o1, o2, chrom, pos):
        '''
        Given a Bedpe object generate the correct alt string for a BND VCF entry
        '''
        alt_string = '{3}{0}{1}:{2}{0}{4}'
        return alt_string.format(self.determine_sep(o2), chrom, pos, *self.determine_flanks(o1))

    def convert(self, bedpe):
        '''
        Convert a bedpe object to Vcf object(s). Returns a list of entries.
        '''
        adjust_tag1, adjust_tag2 = 'CIPOS', 'CIEND'
        if bedpe.malformedFlag == 1:
            adjust_tag2, adjust_tag1 = adjust_tag1, adjust_tag2

        b1 = self.adjust_by_tag(bedpe, adjust_tag1, bedpe.o1, bedpe.s1)
        primary_bedpe_list = [
                bedpe.c1, 
                b1,
                bedpe.orig_name1,
                bedpe.orig_ref1,
                bedpe.orig_alt1,
                bedpe.score,
                bedpe.filter,
                bedpe.info1
                ] + bedpe.misc
        var = Variant(primary_bedpe_list, self.vcf_header)
        to_return = [var]

        if bedpe.svtype == 'BND':
            b2 = self.adjust_by_tag(bedpe, adjust_tag2, bedpe.o2, bedpe.s2)
            
            secondary_bedpe_list = [
                    bedpe.c2,
                    b2,
                    bedpe.orig_name2,
                    bedpe.orig_ref2,
                    bedpe.orig_alt2,
                    bedpe.score,
                    bedpe.filter,
                    bedpe.info2
                    ] + bedpe.misc

            var2 = Variant(secondary_bedpe_list, self.vcf_header)
            if bedpe.malformedFlag == 0:
                to_return += [var2]
            elif bedpe.malformedFlag == 1:
                #Only returning one of our entries
                to_return[0] = var2

        return to_return

