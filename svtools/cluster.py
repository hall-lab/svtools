class Cluster(object):
    '''
    Stores information about overlapping SV and tracks the "best" one of the group
    '''
    def __init__(self):
        self.elements = [None]
        self.chrom_a = None
        self.min_a = float("inf")
        self.max_a = 0
        self.chrom_b = None
        self.min_b = float("inf")
        self.max_b = 0
        self.size = 0
        self.strand_a = ''
        self.strand_b = ''
        self.sv_event = ''
        self.filter = 0

    def can_add(self, bedpe, max_distance):
        '''
        Check whether a bedpe object is addable to this cluster given the max_distance
        '''
        if self.size == 0:
            return True
        if self.sv_event ==  bedpe.svtype:
            if (self.strand_a != bedpe.o1 or
                    self.strand_b != bedpe.o2):
                return False

            if (self.chrom_a != bedpe.c1 or
                    self.min_a - max_distance > bedpe.e1 or
                    self.max_a + max_distance < bedpe.s1):
                return False

            if (self.chrom_b != bedpe.c2 or
                    self.min_b - max_distance > bedpe.e2 or
                    self.max_b + max_distance < bedpe.s2):
                return False
            else:
                return True
        else:
            return False

    def add(self, bedpe, eval_param):
        '''
        Add a new Bedpe object to this cluster
        '''
        if eval_param is None or eval_param.lower() == 'af':
            if  bedpe.af != '.' and bedpe.af > self.filter:
                #First node represents best variant
                self.elements[0] = bedpe
        self.size += 1
        self.sv_event=bedpe.svtype
        self.filter = max(self.filter,bedpe.af)
        self.chrom_a = bedpe.c1
        self.min_a = min(self.min_a, bedpe.s1)
        self.max_a = max(self.max_a, bedpe.e1)
        self.chrom_b = bedpe.c2
        self.min_b = min(self.min_b, bedpe.s2)
        self.max_b = max(self.max_b, bedpe.e2)
        self.strand_a = bedpe.o1
        self.strand_b = bedpe.o2

    def get_cluster_string(self):
        '''
        The string for the best bedpe object for this cluster
        '''
        if self.elements[0] is not None:
            return str(self.elements[0])
        else:
            raise ValueError("Cluster doesn't contain a valid BEDPE line")
