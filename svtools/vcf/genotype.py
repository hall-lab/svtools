import sys

class Genotype(object):
    '''
    This class stores information about each sample.
    '''
    def __init__(self, variant, value_list):
        '''
        Initialize the class. All instances have a GT field, 
        but that is enforced in the Variant class.
        '''
        self.value_list = value_list
        self.variant = variant

    def __eq__(self, other):
        return self.get_gt_string() == other.get_gt_string()
    
    def set_format(self, field, value):
        '''
        Set information for an individual format field.
        '''
        if field in self.variant.format_set:
            if field in self.variant.format_dict:
                i = self.variant.format_dict[field]
                self._set_value(i, value)
            else:
                num_fields = len(self.variant.format_dict)
                self.variant.format_dict[field] = num_fields
                self._set_value(num_fields, value)
        else:
            sys.stderr.write('\nError: invalid FORMAT field, \"' + field + '\"\n')
            sys.exit(1)

    def _set_value(self, index, value):
        try:
            self.value_list[index] = value
        except IndexError:
            self.value_list.extend(['.'] * (index - len(self.value_list)))
            self.value_list.append(value)

    def get_format(self, field):
        '''
        Get value of particular field key
        '''
        try:
            return self.value_list[self.variant.format_dict[field]]
        except IndexError:
            if field != 'GT':
                return '.'
            else:
                return './.'

    def get_gt_string(self):
        '''
        Convert object back to string. 
        
        If some values are missing (at the end for example) they are printed out as 
        all format fields present in any Genotype instance in the Variant line
        are tracked.
        '''
        g_list = list()
        for f in self.variant.format_list:
            if f.id in self.variant.format_dict:
                value = self.get_format(f.id)
                if type(value) == float:
                    g_list.append('%0.2f' % value)
                else:
                    g_list.append(str(value))
        return ':'.join(g_list)

