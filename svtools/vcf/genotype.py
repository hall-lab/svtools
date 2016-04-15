import sys

class Genotype(object):
    '''
    This class stores information about each sample.
    '''
    def __init__(self, variant, gt):
        '''
        Initialize the class. All instances have a GT field.
        '''
        self.format = dict()
        self.variant = variant
        self.set_format('GT', gt)

    def __eq__(self, other):
        return self.get_gt_string() == other.get_gt_string()
    
    def set_formats(self, fields, values):
        '''
        Set many format fields for this instance.
        Updates format information in the owning Variant class.
        '''
        format_set = self.variant.format_set
        add_to_active = self.variant.active_formats.add
        active_formats = self.variant.active_formats
        format_dict = self.format

        for field, value in zip(fields, values):
            if field in format_set:
                format_dict[field] = value
                if field not in active_formats:
                    add_to_active(field)
            else:
                sys.stderr.write('\nError: invalid FORMAT field, \"' + field + '\"\n')
                sys.exit(1)

    def set_format(self, field, value, update_active=True):
        '''
        Set information for an individual format field.
        '''
        if field in self.variant.format_set:
            self.format[field] = value
            if field not in self.variant.active_formats:
                self.variant.active_formats.add(field)
                self.variant.update_active_format_list()
        else:
            sys.stderr.write('\nError: invalid FORMAT field, \"' + field + '\"\n')
            sys.exit(1)

    def get_format(self, field):
        '''
        Get value of particular field key
        '''
        return self.format[field]

    def get_gt_string(self):
        '''
        Convert object back to string. 
        
        If some values are missing (at the end for example) they are printed out as 
        all format fields present in any Genotype instance in the Variant line
        are tracked.
        '''
        g_list = list()
        for f in self.variant.active_format_list:
            if f in self.format:
                if type(self.format[f]) == float:
                    g_list.append('%0.2f' % self.format[f])
                else:
                    g_list.append(str(self.format[f]))
            else:
                g_list.append('.')
        return ':'.join(g_list)

