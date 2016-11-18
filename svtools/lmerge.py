import svtools.l_bp as l_bp
from svtools.breakpoint import Breakpoint
import svtools.logspace as ls

import sys
import numpy as np
import argparse

def null_format_string(format_string):
    null_list = []
    num_null_fields = len(format_string.split(':'))
    if format_string.startswith('GT:'):
        null_list = ['./.']
        num_null_fields -= 1
    null_list.extend(list('.' * num_null_fields))
    null_string = ':'.join(null_list)
    return null_string

def print_var_line(l, genotypes=None):
    A = l.rstrip().split('\t')

    if A[4] == '<INV>' and ('--:0' in A[7] or '++:0' in A[7]):
        [sv_type,chr_l,chr_r,strands,start_l,end_l,start_r,end_r,m] = \
                l_bp.split_v(l)

        STRAND_DICT = dict(x.split(':') for x in m['STRANDS'].split(','))
        for o in STRAND_DICT.keys():
            if STRAND_DICT[o] == '0':
                del(STRAND_DICT[o])
        STRANDS = ','.join(['%s:%s' % (o,STRAND_DICT[o]) for o in STRAND_DICT])

        if STRANDS[:2] == '++':
            ALT = 'N]' + chr_l + ':' + m['END'] + ']'
        elif STRANDS[:2] == '--':
            ALT = '[' + chr_l + ':' + m['END'] + '[N'

        SVTYPE = 'BND'
        CIPOS = m['CIEND']
        CIEND = m['CIPOS']
        CIPOS95 = m['CIEND95']
        CIEND95 = m['CIPOS95']
        IMPRECISE = 'IMPRECISE'
        SU = m['SU']
        PE = m['PE']
        SR = m['SR']
        PRPOS = m['PREND']
        PREND = m['PRPOS']
        SNAME = m['SNAME']
        ALG = m['ALG']
        EVENT = A[2]

        A[4] = ALT
        A[7] = ';'.join(['SVTYPE='   + str(SVTYPE),
                         'STRANDS='  + str(STRANDS),
                         'CIPOS='    + str(CIPOS),
                         'CIEND='    + str(CIEND),
                         'CIPOS95='  + str(CIPOS95),
                         'CIEND95='  + str(CIEND95),
                                       str(IMPRECISE),
                         'SU='       + str(SU),
                         'PE='       + str(PE),
                         'SR='       + str(SR),
                         'PRPOS='    + str(PRPOS),
                         'PREND='    + str(PREND),
                         'ALG='      + str(ALG),
                         'SNAME='    + str(SNAME),
                         'EVENT='    + str(EVENT)])

        # reconstruct the line
        l = '\t'.join(A)

    if A[4] not in ['<DEL>', '<DUP>', '<INV>']:
        [sv_type,chr_l,chr_r,strands,start_l,end_l,start_r,end_r,m] = \
                l_bp.split_v(l)

        CHROM = chr_r
        POS = m['END']
        ID = A[2] + '_2'
        REF = 'N'
        ALT = ''

        if A[4][0] == '[':
            ALT = '[' + chr_l + ':' + A[1] + '[N'
        elif A[4][0] == ']':
            ALT = 'N[' + chr_l + ':' + A[1] + '['
        elif A[4][-1] == '[':
            ALT = ']' + chr_l + ':' + A[1] + ']N'
        elif A[4][-1] == ']':
            ALT = 'N]' + chr_l + ':' + A[1] + ']'

        QUAL = A[5]
        FILTER = '.'
        SVTYPE = 'BND'
        STRANDS = m['STRANDS']
        CIPOS = m['CIEND']
        CIEND = m['CIPOS']
        CIPOS95 = m['CIEND95']
        CIEND95 = m['CIPOS95']
        IMPRECISE = 'IMPRECISE'
        SU = m['SU']
        PE = m['PE']
        SR = m['SR']
        PRPOS = m['PREND']
        PREND = m['PRPOS']
        SNAME = m['SNAME']
        EVENT = A[2]
        ALG = m['ALG']
        SECONDARY = 'SECONDARY'
        MATEID=A[2] + '_1'

        INFO = ';'.join(['SVTYPE='   + str(SVTYPE),
                         'STRANDS='  + str(STRANDS),
                         'CIPOS='    + str(CIPOS),
                         'CIEND='    + str(CIEND),
                         'CIPOS95='  + str(CIPOS95),
                         'CIEND95='  + str(CIEND95),
                                       str(IMPRECISE),
                                       str(SECONDARY),
                         'SU='       + str(SU),
                         'PE='       + str(PE),
                         'SR='       + str(SR),
                         'PRPOS='    + str(PRPOS),
                         'PREND='    + str(PREND),
                         'ALG='      + str(ALG),
                         'SNAME='    + str(SNAME),
                         'EVENT='    + str(EVENT),
                         'MATEID='   + str(MATEID)])

        O = [CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO]

        A[7] += ';MATEID=' + A[2] + '_2'
        A[2] += '_1'
        var1 = '\t'.join(A[:8])
        var2 = '\t'.join([str(o) for o in O])
        if genotypes != None:
            var1 = '\t'.join([var1, genotypes])
            var2 = '\t'.join([var2, genotypes])

        print var1
        print var2

    else:
        var = '\t'.join(A[:8])
        if genotypes != None:
            var = '\t'.join([var, genotypes])
        print var

def merge(BP, sample_order, v_id, use_product, include_genotypes=False):
    if len(BP) == 1:
        A = BP[0].l.rstrip().split('\t')
        #tack on id to SNAME
        s_start=A[7].find('SNAME=')
        s_end=A[7].find(';',s_start)
        sname = None
        if (s_end > -1):
            sname = A[7][s_start + 6:s_end]
            A[7] = A[7][:s_start] + \
                    A[7][s_start:s_end] + \
                    ':' + A[2] + \
                    A[7][s_end:]
        else:
            sname = A[7][s_start + 6:]
            A[7]+= ':' + A[2]

        # reset the id to be unique in this file
        v_id += 1
        A[2] = str(v_id)

        #clip out old mate id
        s_start=A[7].find('MATEID=')
        s_end=A[7].find(';',s_start)
        if (s_end > -1):
            A[7] = A[7][:s_start] + A[7][s_end+1:]
        elif (s_start > -1):
            A[7] = A[7][:s_start]

        #clip out old event id
        s_start=A[7].find('EVENT=')
        s_end=A[7].find(';', s_start)
        if (s_end > -1):
            A[7] = A[7][:s_start] + A[7][s_end+1:]
        elif (s_start > -1):
            A[7] = A[7][:s_start]

        #add new mate
        A[7]+= ';EVENT=' + A[2]

        #add new alg
        if use_product:
            A[7]+= ';ALG=PROD'
        else:
            A[7] += ';ALG=SUM'

        GTS = None
        if include_genotypes:
            null_string = null_format_string(A[8])
            gt_dict = { sname: A[9] }
            GTS = '\t'.join([A[8]] + [gt_dict.get(x, null_string) for x in sample_order])
        print_var_line('\t'.join(A), GTS)
        return v_id

    #Sweep the set.  Find the largest intersecting set.  Remove it.  Continue.
    import heapq

    BP.sort(key=lambda x: x.start_l)

    BP_i = range(len(BP)) # index set of each node in the graph
    C = []

    while len(BP_i) > 0:
        h_l = [] #heap of left breakpoint end coordinates and node id (index). heapq is a min heap and the end coord is what will be used for the sorting.
        max_c = []
        max_c_len = 0
        for i in BP_i:
            # remove anything in the heap that doesn't intersect with the current breakpoint
            while (len(h_l) > 0) and (h_l[0][0] < BP[i].start_l):
                heapq.heappop(h_l)

            heapq.heappush(h_l, (BP[i].end_l, i)) # add to the heap

            # at this point everything in h_l intersects on the left
            # but we need to take into account what is going on on the right
            h_r = [] # heap with rightmost starts
            h_l_i = [x[1] for x in h_l] # this is all of the node ids on the heap currently
            h_l_i.sort(key=lambda x:BP[x].start_r) # sort them by their right start
            for j in h_l_i:
                # remove anything in the heap that doesn't intersect with the current breakpoint on the right end
                while (len(h_r) > 0) and (h_r[0][0] < BP[j].start_r):
                    heapq.heappop(h_r)

                # add something to the right heap
                heapq.heappush(h_r, (BP[j].end_r, j))

                if max_c_len < len(h_r):
                    # max clique! Register what nodes we have
                    max_c_len = len(h_r)
                    max_c = [y[1] for y in h_r]

        C.append(max_c)
        for c in max_c:
            BP_i.remove(c)

    for c in C:
        L = []
        R = []
        for b_i in c:
            b = BP[b_i]
            L.append([b.start_l,b.end_l,b.p_l])
            R.append([b.start_r,b.end_r,b.p_r])

        [start_R, end_R, a_R] = l_bp.align_intervals(R)
        [start_L, end_L, a_L] = l_bp.align_intervals(L)

        p_L = [0] * len(a_L[0])
        p_R = [0] * len(a_R[0])

        for c_i in range(len(c)):
            for i in range(len(a_L[c_i])):
                p_L[i] += a_L[c_i][i]

            for i in range(len(a_R[c_i])):
                p_R[i] += a_R[c_i][i]

        ALG = 'SUM'
        if use_product:
            pmax_i_L = p_L.index(max(p_L))
            pmax_i_R = p_R.index(max(p_R))

            miss = 0
            for c_i in range(len(c)):
                if (a_L[c_i][pmax_i_L] == 0) or (a_R[c_i][pmax_i_R] == 0):
                    miss += 1
            if miss == 0:
                ALG = "PROD"
                ls_p_L = [ls.get_ls(1)] * len(a_L[0])
                ls_p_R = [ls.get_ls(1)] * len(a_R[0])
                for c_i in range(len(c)):
                    for i in range(len(a_L[c_i])):
                        ls_p_L[i] = ls.ls_multiply(ls_p_L[i], ls.get_ls(a_L[c_i][i]))

                    for i in range(len(a_R[c_i])):
                        ls_p_R[i] = ls.ls_multiply(ls_p_R[i], ls.get_ls(a_R[c_i][i]))

                ls_sum_L = ls.get_ls(0)
                ls_sum_R = ls.get_ls(0)

                for ls_p in ls_p_L:
                    ls_sum_L = ls.ls_add(ls_sum_L, ls_p)

                for ls_p in ls_p_R:
                    ls_sum_R = ls.ls_add(ls_sum_R, ls_p)

                p_L = []
                for ls_p in ls_p_L:
                    p_L.append(ls.get_p(ls.ls_divide(ls_p, ls_sum_L)))

                p_R = []
                for ls_p in ls_p_R:
                    p_R.append(ls.get_p(ls.ls_divide(ls_p, ls_sum_R)))

        sum_L = sum(p_L)
        sum_R = sum(p_R)
        p_L = [x/sum_L for x in p_L]
        p_R = [x/sum_L for x in p_R]

        [clip_start_L, clip_end_L] = l_bp.trim(p_L)
        [clip_start_R, clip_end_R] = l_bp.trim(p_R)

        new_start_L = start_L + clip_start_L
        new_end_L = end_L - clip_end_L

        new_start_R = start_R + clip_start_R
        new_end_R = end_R - clip_end_R

        p_L = p_L[clip_start_L:len(p_L)-clip_end_L]
        p_R = p_R[clip_start_R:len(p_R)-clip_end_R]

        s_p_L = sum(p_L)
        s_p_R = sum(p_R)

        p_L = [x/s_p_L for x in p_L]
        p_R = [x/s_p_R for x in p_R]

        max_i_L = p_L.index(max(p_L))
        max_i_R = p_R.index(max(p_R))

        ninefive_i_L_start = max_i_L
        ninefive_i_L_end = max_i_L
        ninefive_i_L_total = p_L[max_i_L]
        updated = 0
        while (ninefive_i_L_total < 0.95):
            if (ninefive_i_L_start <= 0) and (ninefive_i_L_end >= (len(p_L)-1)):
                break
            ninefive_i_L_start = max(0, ninefive_i_L_start - 1)
            ninefive_i_L_end = min(len(p_L)-1, ninefive_i_L_end +1)
            ninefive_i_L_total = sum(p_L[ninefive_i_L_start:ninefive_i_L_end+1])
        ninefive_i_L_start = ninefive_i_L_start - max_i_L
        ninefive_i_L_end = ninefive_i_L_end - max_i_L

        ninefive_i_R_start = max_i_R
        ninefive_i_R_end = max_i_R
        ninefive_i_R_total = p_R[max_i_R]
        updated = 0
        while (ninefive_i_R_total < 0.95):
            if (ninefive_i_R_start <= 0) and (ninefive_i_R_end >= len(p_R)-1):
                break
            ninefive_i_R_start = max(0, ninefive_i_R_start - 1)
            ninefive_i_R_end = min(len(p_R)-1, ninefive_i_R_end +1)
            ninefive_i_R_total = sum(p_R[ninefive_i_R_start:ninefive_i_R_end+1])
        ninefive_i_R_end = ninefive_i_R_end - max_i_R
        ninefive_i_R_start = ninefive_i_R_start - max_i_R
        CIPOS95=str(ninefive_i_L_start) + ',' + str(ninefive_i_L_end)
        CIEND95=str(ninefive_i_R_start) + ',' + str(ninefive_i_R_end)

        CHROM = BP[c[0]].chr_l
        POS = new_start_L + max_i_L
        v_id += 1
        ID = str(v_id)
        REF = 'N'

        ALT = ''
        if BP[c[0]].sv_type == 'BND':
            if BP[c[0]].strands[:2] == '++':
                ALT = 'N]' + \
                        BP[c[0]].chr_r + \
                        ':' + \
                        str(new_start_R + max_i_R) + \
                        ']'
            elif BP[c[0]].strands[:2] == '-+':
                ALT = ']' + \
                        BP[c[0]].chr_r + \
                        ':' + \
                        str(new_start_R + max_i_R) + \
                        ']N'
            elif BP[c[0]].strands[:2] == '+-':
                ALT = 'N[' + \
                        BP[c[0]].chr_r + \
                        ':' + \
                        str(new_start_R + max_i_R) + \
                        '['
            elif BP[c[0]].strands[:2] == '--':
                ALT = '[' + \
                        BP[c[0]].chr_r + \
                        ':' + \
                        str(new_start_R + max_i_R) + \
                        '[N'

        else:
            ALT = '<' + BP[c[0]].sv_type + '>'
        QUAL = 0.0
        FILTER = '.'
        FORMAT = BP[c[0]].l.split('\t')[8]
        SVTYPE = BP[c[0]].sv_type

        STRANDS = ''
        strand_map = {}
        e_type_map = {}

        SU = 0
        PE = 0
        SR = 0

        s_name_list = []

        format_string = None
        gt_dict = dict()

        for b_i in c:
            A = BP[b_i].l.rstrip().split('\t')
            if A[5].isdigit():
                QUAL += float(A[5])

            m = l_bp.to_map(A[7])

            for strand_entry in m['STRANDS'].split(','):
                s_type,s_count = strand_entry.split(':')
                if s_type not in strand_map:
                    strand_map[s_type] = 0
                strand_map[s_type] += int(s_count)

            SU += int(m['SU'])
            PE += int(m['PE'])
            SR += int(m['SR'])

            s_name_list.append(m['SNAME'] + ':' + A[2])

            if include_genotypes:
                if format_string == None:
                    format_string = A[8]

                if format_string == A[8]:
                    gt_dict[m['SNAME']] = A[9]
                else:
                    longer = A[8]
                    shorter = format_string
                    if len(longer) < len(shorter):
                        longer, shorter = shorter, longer

                    if longer.find(shorter) == 0:
                        format_string = longer
                        gt_dict[m['SNAME']] = A[9]
                    else:
                        sys.stderr.write('Unable to merge and include genotypes when FORMAT fields differ across VCF files\n')
                        sys.stderr.write('Previous: {0} Current: {1}\n'.format(format_string, A[8]))
                        sys.stderr.write('Variant: {0}\n'.format(m['SNAME'] + ':' + A[2]))
                        sys.exit(1)

        SNAME=','.join(s_name_list)

        GTS = None
        if include_genotypes:
            null_string = null_format_string(format_string)
            GTS = '\t'.join([format_string] + [gt_dict.get(x, null_string) for x in sample_order])

        strand_types_counts = []
        for strand in strand_map:
            strand_types_counts.append(strand + ':' + str(strand_map[strand]))
        STRANDS = ','.join(strand_types_counts)

        if SVTYPE=='DEL':
            SVLEN = (new_start_L + max_i_L) - (new_start_R + max_i_R)
        else:
            SVLEN = (new_start_R + max_i_R) - (new_start_L + max_i_L)

        # Don't set SVLEN if we have an interchromosomal event. Doesn't make any sense.
        if BP[c[0]].chr_l != BP[c[0]].chr_r:
            SVLEN = None

        END = new_start_R + max_i_R
        CIPOS=','.join([str(x) for x in [-1*max_i_L, len(p_L) - max_i_L - 1]])
        CIEND=','.join([str(x) for x in [-1*max_i_R, len(p_R) - max_i_R - 1]])
        IMPRECISE='IMPRECISE'
        PRPOS=','.join([str(x) for x in p_L])
        PREND=','.join([str(x) for x in p_R])


        if (int(CIPOS.split(',')[0]) > int(CIPOS95.split(',')[0])) or \
            (int(CIPOS.split(',')[1]) < int(CIPOS95.split(',')[1])) or \
            (int(CIEND.split(',')[0]) > int(CIEND95.split(',')[0])) or \
            (int(CIEND.split(',')[1]) < int(CIEND95.split(',')[1])):
            sys.stderr.write(CIPOS + "\t" + str(CIPOS95) + "\n")
            sys.stderr.write(CIEND + "\t" + str(CIEND95) + "\n")

        I = ['SVTYPE='   + str(SVTYPE),
             'STRANDS='  + str(STRANDS)
            ]
        if SVLEN is not None:
            I += ['SVLEN='    + str(SVLEN)]
        I += ['CIPOS='    + str(CIPOS),
             'CIEND='    + str(CIEND),
             'CIPOS95='  + str(CIPOS95),
             'CIEND95='  + str(CIEND95),
                           str(IMPRECISE),
             'SU='       + str(SU),
             'PE='       + str(PE),
             'SR='       + str(SR),
             'PRPOS='    + str(PRPOS),
             'PREND='    + str(PREND),
             'ALG='      + str(ALG),
             'SNAME='    + str(SNAME)]

        if BP[c[0]].sv_type == 'BND':
            I.append('EVENT=' + str(ID))
        else:
            I.append('END=' + str(END))

        INFO = ';'.join(I)

        QUAL = str(QUAL)

        O = [CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO]

        print_var_line('\t'.join([str(o) for o in O]), GTS)
    return v_id

def r_cluster(BP_l, sample_order, v_id, use_product, include_genotypes=False):
    # need to resort based on the right side, then extract clusters
    BP_l.sort(key=lambda x: x.start_r)
    BP_l.sort(key=lambda x: x.chr_r)

    BP_r = []
    BP_max_end_r = -1
    BP_chr_r = ''

    for b in BP_l:
        if (len(BP_r) == 0) or \
           ((b.start_r <= BP_max_end_r) and \
           (b.chr_r == BP_chr_r)):
            BP_r.append(b)
            BP_max_end_r = max(BP_max_end_r, b.end_r)
            BP_chr_r = b.chr_r
        else:
            v_id = merge(BP_r, sample_order, v_id, use_product, include_genotypes)
            BP_r = [b]
            BP_max_end_r = b.end_r
            BP_chr_r = b.chr_r

    if len(BP_r) > 0:
        v_id = merge(BP_r, sample_order, v_id, use_product, include_genotypes)

    return v_id

def l_cluster_by_line(file_name, percent_slop=0, fixed_slop=0, use_product=False, include_genotypes=False):
    v_id = 0
    vcf_lines = []
    vcf_headers = list()
    infile=open(file_name, 'r')

    header = ''
    samples = ''

    for l in infile:
      if l[0] == '#':
        if l[1] != '#':
          samples = l.rstrip().split('\t')[9:]
        else:
          # ignore fileDate
          if l[:10] == '##fileDate':
            continue
          if l not in vcf_headers:
            vcf_headers.append(l)
      if l[0] != '#':
        break

    sample_order = []
    for header in vcf_headers:
      if header[:8] == '##SAMPLE':
        sample_order.append(header.rstrip()[13:-1])
      print header,

    # Add in the sample names etc
    chromosome_header_line = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    if include_genotypes:
        sample_header = '\t'.join(sample_order)
        chromosome_header_line += '\tFORMAT\t' + sample_header
    print chromosome_header_line

    BP_l = []
    BP_sv_type = ''
    BP_max_end_l = -1
    BP_chr_l = ''

    b = Breakpoint(l_bp.parse_vcf_record(l), percent_slop=percent_slop, fixed_slop=fixed_slop)
    BP_l.append(b)
    BP_max_end_l = max(BP_max_end_l, b.end_l)
    BP_chr_l = b.chr_l
    BP_sv_type = b.sv_type

    for l in infile:
      b = Breakpoint(l_bp.parse_vcf_record(l), percent_slop=percent_slop, fixed_slop=fixed_slop)
      if (len(BP_l) == 0) or ((b.start_l <= BP_max_end_l) and (b.chr_l == BP_chr_l) and (b.sv_type == BP_sv_type)):
        BP_l.append(b)
        BP_max_end_l = max(BP_max_end_l, b.end_l)
        BP_chr_l = b.chr_l
        BP_sv_type = b.sv_type
      else:
        v_id = r_cluster(BP_l, sample_order, v_id, use_product, include_genotypes)
        BP_l = [b]
        BP_max_end_l = b.end_l
        BP_sv_type = b.sv_type
        BP_chr_l = b.chr_l

    if len(BP_l) > 0:
        v_id = r_cluster(BP_l, sample_order, v_id, use_product, include_genotypes)

    infile.close()

def description():
    return 'merge LUMPY calls inside a single file from svtools lsort'

def epilog():
    return 'Note that if both slop parameters are set then the maximum is used.'

def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--inFile', metavar='<FILE>', help='a sorted VCF file generated by svtools lsort. Each INFO field must contain an SNAME tag containing the sample name (e.g. SNAME=SAMPLE_NAME)')
    parser.add_argument('-p', '--percent-slop', metavar='<FLOAT>', type=float, default=0.0, help='increase the the breakpoint confidence interval both up and down stream by a given proportion of the original size')
    parser.add_argument('-f', '--fixed-slop', metavar='<INT>', type=int, default=0, help='increase the the breakpoint confidence interval both up and down stream by a given fixed size')
    parser.add_argument('--sum', dest='use_product', action='store_false', default=True, help='calculate breakpoint PDF and position using sum algorithm instead of product')
    parser.add_argument('-g', dest='include_genotypes', action='store_true', default=False, help='include original genotypes in output. When multiple variants are merged, the last will dictate the genotype field')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description(), epilog=epilog())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    l_cluster_by_line(args.inFile,
            percent_slop=args.percent_slop,
            fixed_slop=args.fixed_slop,
            use_product=args.use_product,
            include_genotypes=args.include_genotypes)

if __name__ == "__main__":
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
