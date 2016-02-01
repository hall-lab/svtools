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
                # FIXME This should be an exception
                sys.stderr.write('\nError: invalid FORMAT field, \"' + field + '\"\n')
                sys.exit(1)

    def set_format(self, field, value):
        if field in self.variant.format_set:
            self.format[field] = value
            if field not in self.variant.active_formats:
                self.variant.active_formats.add(field)
        else:
            # FIXME This should be an exception
            sys.stderr.write('\nError: invalid FORMAT field, \"' + field + '\"\n')
            sys.exit(1)

    def get_format(self, field):
        return self.format[field]

    def get_gt_string(self):
       g_list = list()
       for f in self.variant.format_list:
           if f.id in self.variant.active_formats:
               if f.id in self.format:
                   if type(self.format[f.id]) == float:
                       g_list.append('%0.2f' % self.format[f.id])
                   else:
                       g_list.append(str(self.format[f.id]))
               else:
                   g_list.append('.')
       return ':'.join(g_list)

