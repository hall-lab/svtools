import svtools.l_bp as l_bp
from svtools.breakpoint import Breakpoint
import svtools.logspace as ls
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
from svtools.utils import parse_bnd_alt_string, InputStream
from svtools.exceptions import MissingProbabilitiesException

import sys
import numpy as np
import argparse
import heapq
import re

def null_format_string(format_string):
    null_list = []
    num_null_fields = len(format_string.split(':'))
    if format_string.startswith('GT:'):
        null_list = ['./.']
        num_null_fields -= 1
    null_list.extend(list('.' * num_null_fields))
    null_string = ':'.join(null_list)
    return null_string


def merge_single_bp(BP, sample_order, v_id, use_product, vcf, vcf_out, include_genotypes):

    A = BP[0].l.rstrip().split('\t')
    var = Variant(A,vcf)
    try:
        sname = var.get_info('SNAME')
        var.set_info('SNAME', sname + ':' + var.var_id)
    except KeyError:
        pass
    var.var_id=str(v_id)

    if use_product:
        var.set_info('ALG', 'PROD')
    else:
        var.set_info('ALG', 'SUM')

    GTS = None
    if include_genotypes:
        null_string = null_format_string(A[8])
        gt_dict = { sname: A[9] }
        GTS = '\t'.join([gt_dict.get(x, null_string) for x in sample_order])
        var.gts = None
        var.gts_string = GTS

    return var


def order_cliques(BP, C):

    #Sweep the set.  Find the largest intersecting set.  Remove it.  Continue.

    BP_i = range(len(BP)) # index set of each node in the graph
    while len(BP_i) > 0:
        h_l = [] #heap of left breakpoint end coordinates and node id (index). heapq is a min heap and the end coord is what will be used for the sorting.
        max_c = []
        max_c_len = 0
        for i in BP_i:
            # remove anything in the heap that doesn't intersect with the current breakpoint
            while (len(h_l) > 0) and (h_l[0][0] < BP[i].left.start):
                heapq.heappop(h_l)

            heapq.heappush(h_l, (BP[i].left.end, i)) # add to the heap

            # at this point everything in h_l intersects on the left
            # but we need to take into account what is going on on the right
            h_r = [] # heap with rightmost starts
            h_l_i = [x[1] for x in h_l] # this is all of the node ids on the heap currently
            h_l_i.sort(key=lambda x:BP[x].right.start) # sort them by their right start
            for j in h_l_i:
                # remove anything in the heap that doesn't intersect with the current breakpoint on the right end
                while (len(h_r) > 0) and (h_r[0][0] < BP[j].right.start):
                    heapq.heappop(h_r)

                # add something to the right heap
                heapq.heappush(h_r, (BP[j].right.end, j))

                if max_c_len < len(h_r):
                    # max clique! Register what nodes we have
                    max_c_len = len(h_r)
                    max_c = [y[1] for y in h_r]

        C.append(max_c)
        for c in max_c:
            BP_i.remove(c)



def getCI95( p_L, p_R, max_i_L, max_i_R):

    ninefive_i_L_start = max_i_L
    ninefive_i_L_end = max_i_L
    ninefive_i_L_total = p_L[max_i_L]

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
    return [CIPOS95, CIEND95]



def combine_pdfs(BP, c, use_product, weighting_scheme):

    L = []
    R = []
    for b_i in c:
        b = BP[b_i]
        L.append([b.left.start, b.left.end, b.left.p])
        R.append([b.right.start, b.right.end, b.right.p])

    [start_R, end_R, a_R] = l_bp.align_intervals(R)
    [start_L, end_L, a_L] = l_bp.align_intervals(L)

    p_L = [0] * len(a_L[0])
    p_R = [0] * len(a_R[0])
    wts = [1] * len(c)

    for c_i in range(len(c)):

        if weighting_scheme == 'evidence_wt':

            A = BP[c[c_i]].l.rstrip().split('\t', 10)
            m = l_bp.to_map(A[7])
            wt=int(m['SU'])
            #sys.stderr.write("wt\t0\t"+str(wt)+"\n")
            a_L[c_i]=[wt*ali for ali in a_L[c_i]]
            a_R[c_i]=[wt*ari for ari in a_R[c_i]]

        elif weighting_scheme == 'carrier_wt':

            A = BP[c[c_i]].l.rstrip().split('\t', 10)
            m = l_bp.to_map(A[7])
            wt = 1
            if 'SNAME' in m:
                wt=len(m['SNAME'].split(','))
            a_L[c_i]=[wt*ali for ali in a_L[c_i]]
            a_R[c_i]=[wt*ari for ari in a_R[c_i]]

        for i in range(len(a_L[c_i])):
            #sys.stderr.write("L\t"+str(i)+"\t"+str(c_i)+"\t"+str(a_L[c_i][i])+"\n")
            p_L[i] += a_L[c_i][i]

        for i in range(len(a_R[c_i])):
            #sys.stderr.write("R\t"+str(i)+"\t"+str(c_i)+"\t"+str(a_R[c_i][i])+"\n")
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

    [ new_start_L, new_end_L ] = [ start_L + clip_start_L,  end_L - clip_end_L ]
    [ new_start_R, new_end_R ] = [ start_R + clip_start_R, end_R - clip_end_R ]

    p_L = p_L[clip_start_L:len(p_L)-clip_end_L]
    p_R = p_R[clip_start_R:len(p_R)-clip_end_R]

    s_p_L = sum(p_L)
    s_p_R = sum(p_R)

    p_L = [x/s_p_L for x in p_L]
    p_R = [x/s_p_R for x in p_R]

    #sys.exit(1)
    return new_start_L, new_start_R, p_L, p_R, ALG

def create_merged_variant(BP, c, v_id, vcf, use_product, weighting_scheme='unweighted'):

    new_start_L, new_start_R, p_L , p_R, ALG = combine_pdfs(BP, c, use_product, weighting_scheme)


    max_i_L = p_L.index(max(p_L))
    max_i_R = p_R.index(max(p_R))

    [cipos95, ciend95]=getCI95( p_L, p_R, max_i_L, max_i_R)
    new_pos_L = new_start_L + max_i_L
    new_pos_R = new_start_R + max_i_R
    BP0=BP[c[0]]
    A=BP0.l.rstrip().split('\t', 10)

    ALT = ''
    if BP0.sv_type == 'BND':
        if BP0.strands[:2] == '++':
            ALT = 'N]' + BP0.right.chrom + ':' + str(new_pos_R) + ']'
        elif BP0.strands[:2] == '-+':
            ALT =  ']' + BP0.right.chrom + ':' + str(new_pos_R) + ']N'
        elif BP0.strands[:2] == '+-':
            ALT = 'N[' + BP0.right.chrom + ':' + str(new_pos_R) + '['
        elif BP0.strands[:2] == '--':
            ALT =  '[' + BP0.right.chrom + ':' + str(new_pos_R) + '[N'
    else:
        ALT = '<' + BP0.sv_type + '>'

    var_list=[ BP0.left.chrom,
               new_pos_L,
               str(v_id),
               'N',
               ALT,
               0.0,
               '.',
               ''] + A[8:]

    var=Variant(var_list, vcf)

    var.set_info('SVTYPE', BP0.sv_type)
    var.set_info('ALG', ALG)

    if var.get_info('SVTYPE')=='DEL':
        var.set_info('SVLEN', new_pos_L - new_pos_R)
    elif BP0.left.chrom == BP0.right.chrom:
        var.set_info('SVLEN', new_pos_R - new_pos_L)
    else:
        SVLEN = None

    if var.get_info('SVTYPE') == 'BND':
        var.set_info('EVENT', str(v_id))
    else:
        var.set_info('END', new_pos_R )

    var.set_info('CIPOS95', cipos95)
    var.set_info('CIEND95', ciend95)
    var.set_info('CIPOS', ','.join([str(x) for x in [-1*max_i_L, len(p_L) - max_i_L - 1]]))
    var.set_info('CIEND', ','.join([str(x) for x in [-1*max_i_R, len(p_R) - max_i_R - 1]]))
    var.set_info('PRPOS', ','.join([str(x) for x in p_L]))
    var.set_info('PREND', ','.join([str(x) for x in p_R]))

    return var


def combine_var_support(var, BP, c, include_genotypes, sample_order):

    strand_map = {}
    qual = 0.0

    [ SU, PE, SR ] = [0,0,0]

    s_name_list = []
    s1_name_list = []

    format_string = var.get_format_string()
    gt_dict = dict()

    for b_i in c:
        A = BP[b_i].l.rstrip().split('\t')
        if A[5].isdigit():
            qual += float(A[5])

        m = l_bp.to_map(A[7])

        for strand_entry in m['STRANDS'].split(','):
            s_type,s_count = strand_entry.split(':')
            if s_type not in strand_map:
                strand_map[s_type] = 0
            strand_map[s_type] += int(s_count)

        SU += int(m['SU'])
        PE += int(m['PE'])
        SR += int(m['SR'])

        if 'SNAME' in m:
            s_name_list.append(m['SNAME'] + ':' + A[2])

        if include_genotypes:

            if format_string == A[8]:
                gt_dict[m['SNAME']] = A[9]
            else:
                format_dict = dict(zip(A[8].split(':'), A[9].split(':')))
                geno = ':'.join([format_dict.get(i, '.')  for i in var.format_list])
                gt_dict[m['SNAME']] = geno
        else:
            var.format_dict=None

    if s_name_list:
        var.set_info('SNAME', ','.join(s_name_list))

    GTS = None
    if include_genotypes:
        null_string = null_format_string(format_string)
        GTS = '\t'.join([gt_dict.get(x, null_string) for x in sample_order])
        var.gts=None
        var.gts_string=GTS

    strand_types_counts = []
    for strand in strand_map:
        strand_types_counts.append(strand + ':' + str(strand_map[strand]))

    var.set_info('STRANDS', ','.join(strand_types_counts))

    var.qual = qual
    var.set_info('PE', str(PE))
    var.set_info('SU', str(SU))
    var.set_info('SR', str(SR))


def invtobnd(var):

    strands=var.get_info('STRANDS')
    strand_dict = dict(x.split(':') for x in strands.split(','))

    for o in strand_dict.keys():
        if strand_dict[o] == '0':
            del(strand_dict[o])

    strands=','.join(['%s:%s' % (o,strand_dict[o]) for o in strand_dict])
    var.set_info('STRANDS', strands)

    if strands[:2] == '++':
        ALT = 'N]' + var.chrom + ':' + str(var.get_info('END')) + ']'
    elif strands[:2] == '--':
        ALT = '[' + var.chrom + ':' + str(var.get_info('END')) + '[N'

    var.set_info('SVTYPE', 'BND')
    var.alt = ALT

    [ tempci, temp95 ] = [var.get_info('CIPOS'), var.get_info('CIPOS95')]
    try:
        temppr = var.get_info('PRPOS')
    except KeyError:
        raise MissingProbabilitiesException('Required tag PRPOS not found.')

    var.set_info('CIPOS', var.get_info('CIEND'))
    var.set_info('CIEND', tempci)
    var.set_info('CIPOS95', var.get_info('CIEND95'))
    var.set_info('CIEND95', temp95 )
    try:
        var.set_info('PRPOS', var.get_info('PREND'))
    except KeyError:
        raise MissingProbabilitiesException('Required tag PREND not found.')
    var.set_info('PREND', temppr )


def write_var(var, vcf_out, include_genotypes=False):

    v_id=var.var_id
    if var.get_info('CIPOS95') != '0,0' or var.get_info('CIEND95') != '0,0':
        var.set_info('IMPRECISE', True)
    else:
        var.set_info('IMPRECISE', False)

    if var.get_info('SVTYPE') == 'INV' and ('--:0' in var.get_info('STRANDS') or '++:0' in var.get_info('STRANDS')):

        invtobnd(var)

    if var.alt not in ['<DEL>', '<DUP>', '<INV>']:

        var.var_id=str(v_id)+'_1'
        var.set_info('EVENT', v_id)
        var.set_info('MATEID', str(v_id)+'_2')
        var.info.pop('END', None)
        var.info.pop('SVLEN', None)

        varstring=var.get_var_string(use_cached_gt_string=True)
        if not include_genotypes:
            varstring='\t'.join(varstring.split('\t', 10)[:8])

        vcf_out.write(varstring+'\n')

        new_alt = ''

        if var.alt[0] == '[':
            new_alt = '[' + var.chrom + ':' + str(var.pos) + '[N'
        elif var.alt[0] == ']':
            new_alt = 'N[' + var.chrom + ':' + str(var.pos) + '['
        elif var.alt[-1] == '[':
            new_alt = ']' + var.chrom + ':' + str(var.pos) + ']N'
        elif var.alt[-1] == ']':
            new_alt = 'N]' + var.chrom + ':' + str(var.pos) + ']'

        sep, chrom, pos = parse_bnd_alt_string(var.alt)
        var.chrom = chrom
        var.pos = int(pos)
        var.var_id = str(v_id)+'_2'
        var.set_info('MATEID', str(v_id)+'_1')
        var.set_info('SECONDARY', True)
        var.alt = new_alt

        [ tempci, temp95 ] = [var.get_info('CIPOS'), var.get_info('CIPOS95')]
        try:
            temppr = var.get_info('PRPOS')
        except KeyError:
            raise MissingProbabilitiesException('Required tag PRPOS not found.')
        var.set_info('CIPOS', var.get_info('CIEND'))
        var.set_info('CIEND', tempci)
        var.set_info('CIPOS95', var.get_info('CIEND95'))
        var.set_info('CIEND95', temp95 )
        try:
            var.set_info('PRPOS', var.get_info('PREND'))
        except KeyError:
            raise MissingProbabilitiesException('Required tag PREND not found.')
        var.set_info('PREND', temppr )

        varstring=var.get_var_string(use_cached_gt_string=True)
        if not include_genotypes:
            varstring='\t'.join(varstring.split('\t', 10)[:8])

        vcf_out.write(varstring+'\n')


    else:
        varstring=var.get_var_string(use_cached_gt_string=True)
        if not include_genotypes:
            varstring='\t'.join(varstring.split('\t', 10)[:8])

        vcf_out.write(varstring+'\n')



def merge(BP, sample_order, v_id, use_product, vcf, vcf_out, include_genotypes=False, weighting_scheme='unweighted'):

    if len(BP) == 1:
       #merge a single breakpoint
        v_id+=1
        var=merge_single_bp(BP, sample_order, v_id, use_product, vcf, vcf_out, include_genotypes)
        write_var(var, vcf_out, include_genotypes)

    else:

        BP.sort(key=lambda x: x.left.start)
        ordered_cliques = []
        order_cliques(BP, ordered_cliques)

        #merge cliques
        for cliq in ordered_cliques:
            v_id+=1
            var=create_merged_variant(BP, cliq, v_id, vcf, use_product, weighting_scheme)
            combine_var_support(var, BP, cliq, include_genotypes, sample_order)
            write_var(var, vcf_out, include_genotypes)

    return v_id


def r_cluster(BP_l, sample_order, v_id, use_product, vcf, vcf_out, include_genotypes=False, weighting_scheme='unweighted'):

    # need to resort based on the right side, then extract clusters
    BP_l.sort(key=lambda x: x.right.start)
    BP_l.sort(key=lambda x: x.right.chrom)

    BP_r = []
    BP_max_end_r = -1
    BP_chr_r = ''

    for b in BP_l:
        if (len(BP_r) == 0) or \
           ((b.right.start <= BP_max_end_r) and \
           (b.right.chrom == BP_chr_r)):
            BP_r.append(b)
            BP_max_end_r = max(BP_max_end_r, b.right.end)
            BP_chr_r = b.right.chrom
        else:
            v_id = merge(BP_r, sample_order, v_id, use_product, vcf, vcf_out, include_genotypes, weighting_scheme)
            BP_r = [b]
            BP_max_end_r = b.right.end
            BP_chr_r = b.right.chrom

    if len(BP_r) > 0:
        v_id = merge(BP_r, sample_order, v_id, use_product, vcf, vcf_out, include_genotypes, weighting_scheme)

    return v_id



def l_cluster_by_line(file_name, percent_slop=0, fixed_slop=0, use_product=False, include_genotypes=False, weighting_scheme='unweighted'):

    v_id = 0

    in_header = True
    header = []
    vcf = Vcf()
    vcf_out=sys.stdout

    with InputStream(file_name) as vcf_stream:

        BP_l = []
        BP_sv_type = ''
        BP_max_end_l = -1
        BP_chr_l = ''
        sample_order = []

        for line in vcf_stream:

            if in_header:

                if line.startswith('##'):
                    header.append(line)
                    continue

                elif line.startswith('#CHROM'):
                    v=line.rstrip().split('\t')
                    for headline in header:
                        if headline[:8] == '##SAMPLE':
                            sample_order.append(headline.rstrip()[13:-1])
                    hline=''
                    if include_genotypes :
                        v.extend(sample_order)
                        hline='\t'.join(v)
                    else :
                        v=v[:8]
                        hline='\t'.join(v)
                    header.append(hline)
                    in_header=False
                    vcf.add_header(header)
                    vcf.add_info('ALG', '1', 'String', 'Algorithm used to merge this breakpoint')

                    if include_genotypes:
                        vcf_out.write(vcf.get_header()+'\n')
                    else:
                        vcf_out.write(vcf.get_header(False)+'\n')

                continue

            b = Breakpoint(l_bp.parse_vcf_record(line), percent_slop=percent_slop, fixed_slop=fixed_slop)
            if (len(BP_l) == 0) or ((b.left.start <= BP_max_end_l) and (b.left.chrom == BP_chr_l) and (b.sv_type == BP_sv_type)):
                BP_l.append(b)
                BP_max_end_l = max(BP_max_end_l, b.left.end)
                BP_chr_l = b.left.chrom
                BP_sv_type = b.sv_type

            else:
                v_id = r_cluster(BP_l, sample_order, v_id, use_product, vcf, vcf_out, include_genotypes, weighting_scheme)
                BP_l = [b]
                BP_max_end_l = b.left.end
                BP_sv_type = b.sv_type
                BP_chr_l = b.left.chrom

        if len(BP_l) > 0:
            v_id = r_cluster(BP_l, sample_order, v_id, use_product, vcf, vcf_out, include_genotypes, weighting_scheme)

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
    parser.add_argument('-w', dest='weighting_scheme', metavar='<STRING>', default="unweighted", choices=['carrier_wt', 'evidence_wt'], help='weighting scheme (intended for use in tiered merging), options: unweighted, carrier_wt, evidence_wt')
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
            include_genotypes=args.include_genotypes,
            weighting_scheme=args.weighting_scheme)


if __name__ == "__main__":
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
