import sys

class Genotype(object):
    '''
    This class stores information about each sample.
    '''
    def __init__(self, variant, gt):
        '''
        Initialize the class. All instances have a GT field.
        '''
        self.format_dict = dict()
        self.value_list = list()
        self.variant = variant
        self.set_format('GT', gt)

    def __eq__(self, other):
        return self.get_gt_string() == other.get_gt_string()
    
    def initialize_formats(self, field_dict, values):
        '''
        Add formats from field_dict to the object
        '''
        tmp_fields = self.format_dict
        tmp_values = self.value_list
        self.format_dict = field_dict
        self.value_list = values
        for key in tmp_fields:
            if key not in self.format_dict:
                self.format_dict[key] = len(self.value_list)
                self.value_list.append(value)

    def set_format(self, field, value, update_active=True):
        '''
        Set information for an individual format field.
        '''
        if field in self.variant.format_set:
            if field in self.format_dict:
                self.value_list[self.format_dict[field]] = value
            else:
                self.format_dict[field] = len(self.value_list)
                self.value_list.append(value)
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
        return self.value_list[self.format_dict[field]]

    def get_gt_string(self):
        '''
        Convert object back to string. 
        
        If some values are missing (at the end for example) they are printed out as 
        all format fields present in any Genotype instance in the Variant line
        are tracked.
        '''
        g_list = list()
        for f in self.variant.active_format_list:
            if f in self.format_dict:
                value = self.get_format(f)
                if type(value) == float:
                    g_list.append('%0.2f' % value)
                else:
                    g_list.append(str(value))
            else:
                g_list.append('.')
        return ':'.join(g_list)

