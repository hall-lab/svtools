from svtools.vcf.genotype import Genotype
import sys

class Variant(object):
    '''
    Class for storing information stored in a VCF line
    '''
    def __init__(self, var_list, vcf, fixed_genotypes=False):
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
        self.active_formats = set()
        self.active_format_list = list()
        self.gts = dict()
        self.gts_string = None

        # fill in empty sample genotypes
        if len(var_list) < 8:
            sys.stderr.write('\nError: VCF file must have at least 8 columns\n')
            exit(1)

        # make a genotype for each sample at variant
        format_field_tags = var_list[8].split(':')
        self.gts = self._parse_genotypes(format_field_tags, var_list[9:])
        self.update_active_format_list()
        if fixed_genotypes == True:
            self.gts_string='\t'.join(var_list[9:])

        self.info = dict()
        i_split = [a.split('=') for a in var_list[7].split(';')] # temp list of split info column
        for i in i_split:
            if len(i) == 1:
                i.append(True)
            self.info[i[0]] = i[1]

    def _parse_genotypes(self, format_field_tags, genotype_array):
        '''
        Parse the genotype strings
        '''
        gts = dict()
        for index, sample_string in enumerate(genotype_array):
            sample_name = self.sample_list[index]
            try:
                sample_field = sample_string.split(':')
                # sample_name HAS to match the same order.
                gts[sample_name] = Genotype(self, sample_field[0])
                # import the existing fmt fields
                gts[sample_name].set_formats(format_field_tags, sample_field)
            except IndexError:
                gts[sample_name] = Genotype(self, './.')
        return gts

    def update_active_format_list(self):
        '''
        Update the set of this lines 'active' formats.
        This tracks which of the listed formats are actually being used.
        '''
        new_list = list()
        for format in self.format_list:
            if format.id in self.active_formats:
                new_list.append(format.id)
        self.active_format_list = new_list

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

    def get_format_string(self):
        '''
        Construct the FORMAT field containing the names of the fields in the Genotype columns
        '''
        f_list = list()
        for f in self.format_list:
            if f.id in self.active_formats:
                f_list.append(f.id)
        return ':'.join(f_list)

    def get_gt_string(self):
        '''
        Construct the genotype string. If fixed_genotypes was passed in, the cached string is used.
        Otherwise it is constructed from the Genotype objects.
        '''
        if self.gts_string:
            return self.gts_string
        else:
            return '\t'.join(self.genotype(s).get_gt_string() for s in self.sample_list)

    def genotype(self, sample_name):
        '''
        Return the Genotype object for the request sample
        '''
        try:
            return self.gts[sample_name]
        except KeyError as e:
            sys.stderr.write('\nError: invalid sample name, \"' + sample_name + '\"\n')
            sys.exit(1)

    def get_var_string(self, use_cached_gt_string=False):
        '''
        Return the String representation for this line
        '''
        if len(self.active_formats) == 0:
            s = '\t'.join(map(str,[
                self.chrom,
                self.pos,
                self.var_id,
                self.ref,
                self.alt,
                self.qual,
                self.filter,
                self.get_info_string()
            ]))
        elif use_cached_gt_string == False:
            s = '\t'.join(map(str,[
                self.chrom,
                self.pos,
                self.var_id,
                self.ref,
                self.alt,
                self.qual,
                self.filter,
                self.get_info_string(),
                self.get_format_string(),
                '\t'.join(self.genotype(s).get_gt_string() for s in self.sample_list)
            ]))
        else:
            if self.gts_string == None :
                sys.stderr.write("Error no gt_string\n")
                sys.exit(1);
            else:
                s = '\t'.join(map(str,[
                    self.chrom,
                    self.pos,
                    self.var_id,
                    self.ref,
                    self.alt,
                    self.qual,
                    self.filter,
                    self.get_info_string(),
                    self.get_format_string(),
                    self.gts_string
            ]))
        return s
