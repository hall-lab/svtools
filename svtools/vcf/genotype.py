class Genotype(object):
    def __init__(self, variant, sample_name, gt):
        self.format = dict()
        self.variant = variant
        self.set_format('GT', gt)

    def set_format(self, field, value):
        if field in [i.id for i in self.variant.format_list]:
            self.format[field] = value
            if field not in self.variant.active_formats:
                self.variant.active_formats.add(field)
        else:
            # FIXME This should be an exception
            sys.stderr.write('\nError: invalid FORMAT field, \"' + field + '\"\n')
            exit(1)

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
                       g_list.append(self.format[f.id])
               else:
                   g_list.append('.')
       return ':'.join(map(str,g_list))

