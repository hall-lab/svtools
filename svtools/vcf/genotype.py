import sys

class Genotype(object):
    def __init__(self, variant, gt):
        self.format = dict()
        self.variant = variant
        self.set_format('GT', gt)
    
    def set_formats(self, fields, values):
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
        if field in self.variant.format_set:
            self.format[field] = value
            if field not in self.variant.active_formats:
                self.variant.active_formats.add(field)
                self.variant.update_active_format_list()
        else:
            sys.stderr.write('\nError: invalid FORMAT field, \"' + field + '\"\n')
            sys.exit(1)

    def get_format(self, field):
        return self.format[field]

    def get_gt_string(self):
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

