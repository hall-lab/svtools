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
    def adjust_by_cipos(bedpe):
        '''
        Undo adjustments to BEDPE coordinates for the VCF POS based on confidence intervals
        '''
        position = bedpe.s1
        if bedpe.o1 == '-' and bedpe.svtype == 'BND':
            position += 1   #undo left adjust based on strandedness
        if 'CIPOS=' in bedpe.info1:
            cipos = re.split('=|,', ''.join(filter(lambda x: 'CIPOS=' in x, bedpe.info1.split(';'))))
            position -= int(cipos[1])
        return position

    @staticmethod
    def adjust_by_ciend(bedpe):
        '''
        Undo adjustments to BEDPE coordinates for the VCF END field based on confidence intervals
        '''
        end = bedpe.s2
        if bedpe.o2 == '-' and bedpe.svtype == 'BND':
            end += 1
    
        if 'CIEND=' in bedpe.info1:     
            ciend = re.split('=|,', ''.join(filter(lambda x: 'CIEND=' in x, bedpe.info1.split(';'))))
            end -= int(ciend[1])
        return end

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
        b1 = self.adjust_by_cipos(bedpe)
        primary_bedpe_list = [
                bedpe.c1, 
                b1,
                bedpe.name,
                'N',
                '<' + str(bedpe.svtype) + '>', #ALT
                bedpe.score,
                bedpe.filter,
                bedpe.info1
                ] + bedpe.misc
        var = Variant(primary_bedpe_list, self.vcf_header)
        to_return = [var]

        if bedpe.svtype == 'BND':
            b2 = self.adjust_by_ciend(bedpe)
            var.var_id += '_1'
            var.alt = self.bnd_alt_string(bedpe.o1, bedpe.o2, bedpe.c2, b2)
            
            secondary_bedpe_list = [
                    bedpe.c2,
                    b2,
                    bedpe.name + '_2',
                    'N',
                    self.bnd_alt_string(bedpe.o2, bedpe.o1, bedpe.c1, b1),
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

