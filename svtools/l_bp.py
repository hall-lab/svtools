import re

def find_all(a_str, sub):
    '''
    Generator to find the coordinates of the start of all
    occurrences of a substring
    '''
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches

def parse_vcf(vcf_file_stream, vcf_lines, vcf_headers, add_sname=True, include_ref=False):
    header = ''
    samples = ''

    for l in vcf_file_stream:
        if l[0] == '#':
            if l[1] != '#':
                samples = l.rstrip().split('\t')[9:]
            else:
                # ignore fileDate
                if l[:10] == '##fileDate':
                    continue
                if l not in vcf_headers:
                    vcf_headers.append(l)
        else:
            A = l.split('\t')
            if not include_ref and ('GT' in A[8]):
                has_nonref = False
                for sample_field in A[9:]:
                    if not (sample_field.startswith('0/0') or sample_field.startswith('./.')):
                        has_nonref = True
                        break
                if not has_nonref:
                    continue

            if not 'SECONDARY' in A[7]:

                if add_sname and (samples != ''):
                    A[7] += ';' + 'SNAME=' + ','.join(samples)
                    l = '\t'.join(A)


                if 'SVTYPE=BND' in A[7]:
                    m = re.search(r"(\[|\])(.*)(\[|\])",A[4])
                    o_chr,o_pos = m.group(2).split(':')

                    if (o_chr == A[0]) and (('--:' in A[7]) != ('++' in A[7])):
                        neg_s = A[7].find('--:')
                        pos_s = A[7].find('++:')

                        if neg_s > 0:
                            neg_e = neg_s + A[7][neg_s:].find(';')
                            pre=A[7][:neg_s]
                            mid=A[7][neg_s:neg_e]
                            post=A[7][neg_e:]
                            A[7] = pre + '++:0,' + mid + post
                        else:
                            pos_e = pos_s + A[7][pos_s:].find(';')
                            pre=A[7][:pos_s]
                            mid=A[7][pos_s:pos_e]
                            post=A[7][pos_e:]
                            A[7] = pre + mid + ',--:0' + post

                        A[7] = 'SVTYPE=INV' + A[7][10:] + ';END=' + o_pos
                        A[4] = '<INV>'
                        l = '\t'.join(A)
                vcf_lines.append(l)

    return samples

def parse_vcf_record(vcf_line):

    A = vcf_line.split('\t')
    if not 'SECONDARY' in A[7]:

        if 'SVTYPE=BND' in A[7]:
            m = re.search(r"(\[|\])(.*)(\[|\])",A[4])
            o_chr, o_pos=m.group(2).split(':')

            if (o_chr == A[0]) and (('--:' in A[7]) != ('++' in A[7])):
                neg_s = A[7].find('--:')
                pos_s = A[7].find('++:')

                if neg_s > 0:
                    neg_e = neg_s + A[7][neg_s:].find(';')
                    pre=A[7][:neg_s]
                    mid=A[7][neg_s:neg_e]
                    post=A[7][neg_e:]
                    A[7] = pre + '++:0,' + mid + post
                else:
                    pos_e = pos_s + A[7][pos_s:].find(';')
                    pre=A[7][:pos_s]
                    mid=A[7][pos_s:pos_e]
                    post=A[7][pos_e:]
                    A[7] = pre + mid + ',--:0' + post

                A[7] = 'SVTYPE=INV' + A[7][10:] + ';END=' + o_pos
                A[4] = '<INV>'
                vcf_line='\t'.join(A)

    return vcf_line

def split_v(l):
    '''
    Split a VCF line into constituents and return a subset of values in an array
    '''
    A = l.split('\t', 8)
    m = to_map(A[7])

    chr_l = A[0]
    pos_l = int(A[1])

    chr_r = A[0]
    pos_r = int(A[1])
    if m['SVTYPE'] == 'BND':
        sep = '['
        if not sep in A[4]:
            sep = ']'
        s,e = [x for x in find_all(A[4],sep)]
        chr_r,pos_r = A[4][s+1:e].split(':')
        m['END'] = pos_r
        pos_r = int(pos_r)
    else:
        pos_r = int(m['END'])

    start_l = pos_l + int(m['CIPOS'].split(',')[0])
    end_l = pos_l + int(m['CIPOS'].split(',')[1])

    start_r = pos_r + int(m['CIEND'].split(',')[0])
    end_r = pos_r + int(m['CIEND'].split(',')[1])

    strands = m['STRANDS']

    return [m['SVTYPE'],chr_l,chr_r,strands,start_l,end_l,start_r,end_r,m]

def to_map(s):
    '''
    Convert a string like one would find in a VCF info field to a dictionary
    '''
    m = {}
    for k_v in s.split(';'):
        A = k_v.split('=')
        if len(A) > 1:
            m[A[0]] = A[1]
        else:
            m[A[0]] = None

    return m

def vcf_line_key(l1):
    v1 = split_v(l1)[:8]
    v1[3] = v1[3][:2]
    return v1


def vcf_line_cmp(l1, l2):
    v1 = split_v(l1)
    v2 = split_v(l2)

    v1[3] = v1[3][:2]
    v2[3] = v2[3][:2]

    for i in range(len(v1)-1):
        if v1[i] != v2[i]:
            return cmp(v1[i],v2[i])
    return 0

def header_line_cmp(l1, l2):
    order = ['##source', \
             '##INFO', \
             '##ALT', \
             '##FORMAT',\
             '##SAMPLE']

    # make sure ##fileformat is first
    if l1[:12] == '##fileformat':
        return -1

    if l2[:12] == '##fileformat':
        return 1

    # make sure #CHROM ... is last
    if l1[1] != '#':
        return 1
    elif l2[1] != '#':
        return -1

    if l1.find('=') == -1:
        return -1
    if l2.find('=') == -1:
        return 1

    h1 = l1[:l1.find('=')]
    h2 = l2[:l2.find('=')]
    if h1 not in order:
        return -1
    if h2 not in order:
        return 1
    return cmp(order.index(h1),order.index(h2))

def trim(A):
    '''
    Return offset of first non-zero value from each end of an array
    '''
    clip_start = 0
    for i in range(len(A)):
        if A[i] == 0:
            clip_start += 1
        else:
            break
    clip_end = 0
    for i in range(len(A)-1,-1,-1):
        if A[i] == 0:
            clip_end += 1
        else:
            break
    return [clip_start, clip_end]


# I has 3 components [[start],[end],[p array]]
def align_intervals(I):
    '''
    Find range containing all intervals then pad out
    probabilities with zeroes so each set is of the same length
    '''
    start = -1
    end = -1
    new_I = []

    START = 0
    END = 1
    P = 2

    # find ends
    for i in I:
        if start == -1:
            start = i[START]
            end = i[END]
        else:
            if i[START] < start:
                start = i[START]

            if i[END] > end:
                end = i[END]

    for i in I:
        new_i = i[P]

        if i[START] > start:
            n = i[START] - start
            new_i = [0]*n + new_i

        if i[END] < end:
            n = end - i[END]
            new_i = new_i + [0]*n

        new_I.append(new_i)

    # one interval. last element is array of arrays of probs covering entire interval
    return [start, end, new_I]
