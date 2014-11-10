#!/usr/bin/env python
import sys
import numpy as np
import optparse
from optparse import OptionParser


parser = OptionParser()

parser.add_option("-b",
    "--bedpe_file",
    dest="bedpe_file",
    help="BEDPE file")

parser.add_option("-c",
    "--config_file",
    dest="config_file",
    help="Tab-delim sample config file of NAME id TYPE. " + \
            "Example:NA12878 10  PE")

(options, args) = parser.parse_args()

if not options.bedpe_file:
    parser.error('BEDPE file not given')

if not options.config_file:
    parser.error('config file not given')

C_by_id = {}
C_id_order = {}

f = open(options.config_file,'r')
for l in f:
    A = l.rstrip().split()
    name = A[0]
    sample_id = A[1]
    data_type = A[2]

    C_by_id[sample_id] = [name,data_type]

    C_id_order[sample_id] = len(C_id_order)

f.close()

f = open(options.bedpe_file,'r')

for l in f:
    A = l.rstrip().split('\t')

    who=[]
    name_type_amount = []
    types=[]
    total_amount = 0
    amount_of_each = [0]*len(C_id_order)

    sv_type = A[10][5:8]

    if (A[0] == A[3]) and (int(A[5]) - int(A[1])<1000000):
        sv_type = 'LOCAL_' + sv_type
    else:
        sv_type = 'DISTANT_' + sv_type

    l_width = int(A[2]) - int(A[1])
    r_width = int(A[5]) - int(A[4])

    for sample_id_amount in A[11][4:].split(';'):
        sample_id,amount = sample_id_amount.split(',')

        total_amount += int(amount)
        amount_of_each[C_id_order[sample_id]] = int(amount)

        if not C_by_id[sample_id][0] in who:
             who.append(C_by_id[sample_id][0])

        name_type_amount.append(\
                C_by_id[sample_id][0] + '_' + C_by_id[sample_id][1] + ':' +\
                amount)

        if not C_by_id[sample_id][1] in types:
            types.append(C_by_id[sample_id][1])

    type_string = ''
    if len(types) == 1:
        type_string = types[0] + '_ALONE'
    else:
        type_string = '_'.join(types)

    print '\t'.join( \
            A[0:7] + \
            [\
                str(total_amount), \
                A[8],
                A[9],
                sv_type,
                str(l_width),
                str(r_width),
                A[7],
                ','.join(who), \
                ','.join(name_type_amount),
                A[12],
                type_string,\
                str(len(who)),\
                str(sum(amount_of_each))\
            ] + \
            [str(x) for x in amount_of_each]\
            )

f.close()
