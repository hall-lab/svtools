#!/usr/bin/env python

'''
This script intended to read in an lsort output file and reorders each format field such that CN is the last field, but other ordering is unchanged.
This is useful as a pre-processing step for allowing merging of genotypes in lmerge with v0.3.1
'''

import sys
import itertools

if __name__ == '__main__':
    for line in sys.stdin:
        if line.startswith('#'):
            sys.stdout.write(line)
        else:
            fields = line.rstrip().split('\t')
            format_field = fields[8]
            sample_field = fields[9]
            value_for = dict()
            for key, value in itertools.izip_longest(format_field.split(':'), sample_field.split(':'), fillvalue='.'):
                value_for[key] = value
            new_keys = [ x for x in format_field.split(':') if x != 'CN' and x != 'ASC']
            new_values = [ value_for[x] for x in format_field.split(':') if x != 'CN' and x != 'ASC']
            if 'ASC' in value_for:
                new_keys.append('ASC')
                new_values.append(value_for['ASC'])
            if 'CN' in value_for:
                new_keys.append('CN')
                new_values.append(value_for['CN'])
            sys.stdout.write('\t'.join([ '\t'.join(fields[:8]), ':'.join(new_keys), ':'.join(new_values)]))
            sys.stdout.write('\n')


