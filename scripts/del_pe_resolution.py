#!/usr/bin/env python

from __future__ import division
import json
from collections import Counter
import argparse


def calc_insert_density(hist):
    '''Transform a histogram of counts to a density'''
    dens = Counter()
    total = sum(hist.values())
    for i in list(hist):
        dens[i] = float(hist[i])/total
    return dens


def overlap(dens, shift):
    '''Shift a density over by "shift" bp and calculate the overlap'''
    total = 0.0
    for x in xrange(shift, max(dens) + 1):
        total += min(dens[x], dens[x - shift])
    return total


def find_overlap(dens, target):
    '''Find amount to shift the density to achieve the target overlap value'''
    shift = max(dens) - 1
    current = overlap(dens, shift)
    last = current
    while shift >= 0 and current <= target:
        last = current
        shift -= 1
        current = overlap(dens, shift)
    return (shift + 1, last)


def load_svtyper_json(json_file):
    '''Load an svtyper json'''
    with open(json_file) as f:
        doc = json.load(f)
    return doc


def create_hist(lib):
    '''Create a histogram from svtyper json information'''
    return Counter({
        int(k): int(v) for k, v in lib['histogram'].items()
        })


def calculate_overlaps(doc, target):
    '''Calculate the minimum variant size with target discriminating power'''
    for sample in doc:
        for lib in doc[sample]['libraryArray']:
            hist = create_hist(lib)
            dens = calc_insert_density(hist)
            (size, overlap_prob) = find_overlap(dens, target)
            return (sample, size, overlap_prob)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculate variant size resolution based on cutoff',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('json_file', nargs='+',
                        help='svtyper json files to evaluate'
                        )
    parser.add_argument('--overlap', type=float, metavar='FLOAT',
                        default=0.05,
                        help='maximum density of overlap between discordant \
                        and concordant insert size distributions'
                        )
    args = parser.parse_args()
    print '\t'.join(('Sample', 'MinimumSize', 'Overlap'))
    for f in args.json_file:
        doc = load_svtyper_json(f)
        results = calculate_overlaps(doc, 0.05)
        print '\t'.join([str(x) for x in results])


# Here thar be tests
def test_calc_insert_density():
    t = Counter({1: 1, 2: 2, 3: 1})
    expected = Counter({1: 0.25, 2: 0.5, 3: 0.25})
    assert(calc_insert_density(t) == expected)


def test_overlap():
    t = Counter({1: 0.25, 2: 0.5, 3: 0.25})
    assert(overlap(t, 0) == 1.0)
    assert(overlap(t, 3) == 0.0)
    assert(overlap(t, 1) == 0.5)


def test_find_overlap():
    t = Counter({1: 0.2, 2: 0.5, 3: 0.3})
    assert(find_overlap(t, 1.0) == (0, 1.0))
    assert(find_overlap(t, 0.5) == (1, 0.5))
