from svtools.vcf.genotype import Genotype
import sys

class Variant(object):
    def __init__(self, var_list, vcf, fixed_genotypes=False):
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
            # FIXME This should be an exception
            sys.stderr.write('\nError: VCF file must have at least 8 columns\n')
            exit(1)

        # make a genotype for each sample at variant
        format_field_tags = var_list[8].split(':')
        for s in self.sample_list:
            try:
                sample_field = var_list[vcf.sample_to_col(s)].split(':')
                self.gts[s] = Genotype(self, sample_field[0])
                # import the existing fmt fields
                self.gts[s].set_formats(format_field_tags, sample_field)
            except IndexError:
                self.gts[s] = Genotype(self, './.')
        self.update_active_format_list()
        if fixed_genotypes == True:          
            self.gts_string='\t'.join(var_list[9:])

        self.info = dict()
        i_split = [a.split('=') for a in var_list[7].split(';')] # temp list of split info column
        for i in i_split:
            if len(i) == 1:
                i.append(True)
            self.info[i[0]] = i[1]

    def update_active_format_list(self):
        new_list = list()
        for format in self.format_list:
            if format.id in self.active_formats:
                new_list.append(format.id)
        self.active_format_list = new_list

    def set_info(self, field, value):
        if field in [i.id for i in self.info_list]:
            self.info[field] = value
        else:
            # FIXME This should be an exception
            sys.stderr.write('\nError: invalid INFO field, \"' + field + '\"\n')
            exit(1)

    def get_info(self, field):
        return self.info[field]

    def get_info_string(self):
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
        f_list = list()
        for f in self.format_list:
            if f.id in self.active_formats:
                f_list.append(f.id)
        return ':'.join(f_list)
    
    def get_gt_string(self):
        if self.gts_string:
            return self.gts_string
        else:
            return '\t'.join(self.genotype(s).get_gt_string() for s in self.sample_list)

    def genotype(self, sample_name):
        try:
            return self.gts[sample_name]
        except KeyError as e:
            sys.stderr.write('\nError: invalid sample name, \"' + sample_name + '\"\n')
            sys.exit(1)

    def get_var_string(self, use_cached_gt_string=False):
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
