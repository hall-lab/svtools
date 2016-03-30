import re
import time

class Vcf(object):
    def __init__(self):
        self.file_format = 'VCFv4.2'
        self.other_meta = []
        self.sample_list = []
        self.sample_indices = dict()
        self.info_list = []
        self.format_list = []
        self.filter_list = []
        self.alt_list = []
        self.add_format('GT', 1, 'String', 'Genotype')

    def parse_meta(self, line):
        a = line[line.find('<')+1:line.rfind('>')]
        r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
        return r.findall(a)

    def add_header(self, header):
        for line in header:
            split_header = line.split('=')
            if split_header[0] == '##fileformat':
                self.file_format = line.rstrip().split('=')[1]
            elif split_header[0] == '##fileDate':
                self.other_meta.append('##fileDate=' + time.strftime('%Y%m%d'))
            elif split_header[0] == '##INFO':
                self.add_info(*[b.split('=')[1] for b in self.parse_meta(line)])
            elif split_header[0] == '##ALT':
                self.add_alt(*[b.split('=')[1] for b in self.parse_meta(line)])
            elif split_header[0] == '##FORMAT':
                self.add_format(*[b.split('=')[1] for b in self.parse_meta(line)])
            elif split_header[0] == '##FILTER':
                self.add_filter(*[b.split('=')[1] for b in self.parse_meta(line)])
            elif split_header[0].startswith('##'):
                self.other_meta.append(line.rstrip())
            elif line[0] == '#' and line[1] != '#':
                self.sample_list = line.rstrip().split('\t')[9:]
                for i in xrange(0, len(self.sample_list)):
                    self.sample_indices[self.sample_list[i]] = i + 9

    # return the VCF header
    def get_header(self, include_samples=True):
        header_meta_array = ['##fileformat=' + self.file_format] + \
                        self.other_meta + \
                        [i.hstring for i in self.info_list] + \
                        [a.hstring for a in self.alt_list] + \
                        [f.hstring for f in self.filter_list] + \
                        [f.hstring for f in self.format_list]
        header_array = ['#CHROM',
                        'POS',
                        'ID',
                        'REF',
                        'ALT',
                        'QUAL',
                        'FILTER',
                        'INFO']
        if include_samples:
            header_array += ['FORMAT'] + self.sample_list
        return '\n'.join(header_meta_array + ['\t'.join(header_array)])

    def add_info(self, id, number, type, desc):
        if id not in [i.id for i in self.info_list]:
            inf = self.Info(id, number, type, desc)
            self.info_list.append(inf)

    def add_info_after(self, insert_id, id, number, type, desc):
        for i in range(0, len(self.info_list)):
            if insert_id == self.info_list[i].id and (id not in [j.id for j in self.info_list]):
                inf = self.Info(id, number, type, desc)
                self.info_list.insert(i + 1, inf)
                return

    def add_alt(self, id, desc):
        if id not in [a.id for a in self.alt_list]:
            alt = self.Alt(id, desc)
            self.alt_list.append(alt)

    def add_format(self, id, number, type, desc):
        if id not in [f.id for f in self.format_list]:
            fmt = self.Format(id, number, type, desc)
            self.format_list.append(fmt)

    def add_filter(self, id, desc):
        if id not in [f.id for f in self.filter_list]:
            flt = self.Filter(id, desc)
            self.filter_list.append(flt)

    def add_sample(self, name):
        self.sample_list.append(name)
        # XXX Probably slow. Optimize if we start adding lots of samples ever
        self.sample_indices[name] = self.sample_list.index(name) + 9

    # get the VCF column index of a sample
    # NOTE: this is zero-based, like python arrays
    def sample_to_col(self, sample):
        return self.sample_indices[sample]

    class Info(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##INFO=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'
        def __eq__(self, other):
            return self.hstring == other.hstring

    class Alt(object):
        def __init__(self, id, desc):
            self.id = str(id)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##ALT=<ID=' + self.id + ',Description=\"' + self.desc + '\">'
        def __eq__(self, other):
            return self.hstring == other.hstring

    class Format(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##FORMAT=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'
        def __eq__(self, other):
            return self.hstring == other.hstring

    class Filter(object):
        def __init__(self, id, desc):
            self.id = str(id)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##FILTER=<ID=' + self.id + ',Description=\"' + self.desc + '\">'
        def __eq__(self, other):
            return self.hstring == other.hstring
