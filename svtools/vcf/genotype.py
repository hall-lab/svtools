import sys

class Genotype(object):
    def __init__(self, variant, sample_name, gt):
        self.format = dict()
        self.variant = variant
        self.set_format('GT', gt)

    def set_formats(self, fields, values):
        for field, value in zip(fields, values):
            if field in self.variant.format_set:
                self.format[field] = value
                if field not in self.variant.active_formats:
                    self.variant.active_formats.add(field)
            else:
                sys.stderr.write('\nError: invalid FORMAT field, \"' + field + '\"\n')
                sys.exit(1)

    def set_format(self, field, value):
        self.set_formats( [field], [value])
        self.variant.update_active_format_list()

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

