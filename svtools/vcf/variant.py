from svtools.vcf.genotype import Genotype
import sys

class Variant(object):
    '''
    Class for storing information stored in a VCF line
    '''
    def __init__(self, var_list, vcf):
        '''
        Initialize values.

        If fixed_genotypes is True then a string corresponding to the
        genotype portion of the line is cached for printing later.
        '''
        self.chrom = var_list[0]
        self.pos = int(var_list[1])
        self.var_id = var_list[2]
        self.ref = var_list[3]
        self.alt = var_list[4]
        self.qual = var_list[5]
        self.filter = var_list[6]
        self.sample_list = vcf.sample_list
        self.info_list = vcf.info_list
        self.info = dict()
        self.format_list = vcf.format_list
        self.format_set = {i.id for i in vcf.format_list}
        self.gts = None

        # fill in empty sample genotypes
        if len(var_list) < 8:
            sys.stderr.write('\nError: VCF file must have at least 8 columns\n')
            exit(1)

        # make a genotype for each sample at variant
        self.format_string = var_list[8]
        self.format_dict = { key: index for index, key in enumerate(self.format_string.split(':')) }
        self.gts_string = '\t'.join(var_list[9:])

        if 'GT' not in self.format_dict:
            self.format_dict['GT'] = len(self.format_dict) #add GT if it doesn't exist
            self._uncache_gts()

        self.info = dict()
        i_split = [a.split('=') for a in var_list[7].split(';')] # temp list of split info column
        for i in i_split:
            if len(i) == 1:
                i.append(True)
            self.info[i[0]] = i[1]

    def _parse_genotypes(self, genotype_array):
        '''
        Parse the genotype strings
        '''
        gts = dict()
        for index, sample_string in enumerate(genotype_array):
            sample_name = self.sample_list[index]
            sample_field = sample_string.split(':')
            g = Genotype(self, sample_field)
            gts[sample_name] = g
        return gts

    def set_info(self, field, value):
        '''
        Set value of the specified field in the INFO section.
        The INFO field must exist already.
        '''
        if field in [i.id for i in self.info_list]:
            self.info[field] = value
        else:
            sys.stderr.write('\nError: invalid INFO field, \"' + field + '\"\n')
            exit(1)

    def get_info(self, field):
        '''
        Get a value for the given INFO field
        '''
        return self.info[field]

    def get_info_string(self):
        '''
        Construct the INFO string for printing. Order is matched to the header.
        '''
        i_list = list()
        for info_field in self.info_list:
            if info_field.id in self.info.keys():
                if info_field.type == 'Flag':
                    if self.info[info_field.id]:
                        i_list.append(info_field.id)
                else:
                    i_list.append('%s=%s' % (info_field.id, self.info[info_field.id]))
        return ';'.join(i_list)

    def get_format_string(self, use_cached_format_string=False):
        '''
        Construct the FORMAT field containing the names of the fields in the Genotype columns
        '''
        if use_cached_format_string or self.gts is None:
            return self.format_string
        else:
            f_list = list()
            for f in self.format_list:
                if f.id in self.format_dict:
                    f_list.append(f.id)
            return ':'.join(f_list)

    def get_gt_string(self, use_cached_gt_string=False):
        '''
        Construct the genotype string.
        '''
        if self.gts:
            if use_cached_gt_string:
                return self.gts_string
            else:
                return '\t'.join(self.genotype(s).get_gt_string() for s in self.sample_list)
        else:
            return self.gts_string

    def _uncache_gts(self):
        '''
        Parse genotypes if they are requested
        '''
        if self.gts is None:
            self.gts = self._parse_genotypes(self.gts_string.split('\t'))

    def genotypes(self):
        '''
        Return a list of all genotype data in the Variant line
        '''
        self._uncache_gts()
        return self.gts.values()

    def genotype(self, sample_name):
        '''
        Return the Genotype object for the requested sample
        '''
        self._uncache_gts()
        try:
            return self.gts[sample_name]
        except KeyError:
            sys.stderr.write('\nError: invalid sample name, \"' + sample_name + '\"\n')
            sys.exit(1)

    def set_genotype(self, sample_name, new_genotype):
        '''
        Set the Genotype object for the given sample. Programmer needs to be
        very careful about what gets added here as there is no error checking.
        '''
        self._uncache_gts()
        try:
            self.gts[sample_name] = new_genotype
        except KeyError:
            sys.stderr.write('\nError: invalid sample name, \"' + sample_name + '\"\n')
            sys.exit(1)

    def get_var_string(self, use_cached_gt_string=False):
        '''
        Return the String representation for this line
        '''
        fields = [
                self.chrom,
                self.pos,
                self.var_id,
                self.ref,
                self.alt,
                self.qual,
                self.filter,
                self.get_info_string()
                ]

        if self.format_dict:
            gts_string = self.get_gt_string(use_cached_gt_string)
            if gts_string is None:
                sys.stderr.write("Unable to construct or retrieve genotype string\n")
                sys.exit(1)
            else:
                fields += [
                        self.get_format_string(use_cached_gt_string),
                        gts_string
                    ]
        return '\t'.join(map(str, fields))
