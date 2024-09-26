#!/usr/bin/env python3

import argparse
import sys
import csv
import itertools
import operator
import numpy as np

import common


def common_prefix(s1, s2):
    for i, (x, y) in enumerate(zip(s1, s2)):
        if x != y:
            return i
    return min(len(s1), len(s2))


class Alleles:
    def __init__(self, s):
        self.s = s
        if s.startswith('<'):
            self.s = s
            self.options = ((None,) * 4,)
            return

        if '/' in s:
            pref, suff = s.rsplit(':', 1)
            s = '|'.join(f'{pref}:{subsuff}' for subsuff in suff.split('/'))

        self.options = set()
        for allele in s.split('|'):
            assert ',' not in allele and ';' not in allele and len(allele) > 1,   repr(allele)
            if allele.endswith('*'):
                allele = allele[:-1]

            if ':' not in allele:
                first = int(allele[:3])
                second = int(allele[3:5]) if len(allele) > 3 else None
                third = int(allele[5:].rstrip('NLQ')) if len(allele) > 5 else None
                self.options.add((first, second, third))
            else:
                opt = []
                for field in allele.split(':'):
                    if field[-1] in 'NLQ':
                        field = field[:-1]
                    if field == 'new':
                        break
                    opt.append(int(field))
                self.options.add(tuple(opt))

    def __str__(self):
        return self.s

    def __repr__(self):
        return self.s

    def accuracy(self, oth):
        m = -1
        l = 0
        for t1 in self.options:
            for t2 in oth.options:
                a = common_prefix(t1, t2)
                if a > m:
                    m = a
                    l = len(t1)
        return m, l


def load_alleles(f):
    d = {}
    for line in f:
        if line.startswith('#') or line.startswith('gene\t'):
            continue

        line = line.strip().split('\t')
        gene = line[0]
        sample = line[1]

        alleles = []
        for entry in line[2:]:
            entry = entry.split(' ')
            if len(entry) > 1:
                entry, qual = entry
            else:
                entry = entry[0]
                qual = None
            for allele in entry.split(';'):
                alleles.append((Alleles(allele), qual))
        d[(gene, sample)] = alleles
    return d


def permute_best(alleles1, alleles2):
    if len(alleles1) >= len(alleles2):
        perm_alleles1 = tuple(itertools.permutations(alleles1))
        perm_alleles2 = (alleles2,)
    else:
        perm_alleles1 = (alleles1,)
        perm_alleles2 = tuple(itertools.permutations(alleles2))

    best_acc = -1
    best_i = 0
    best_j = 0
    for i, alleles1 in enumerate(perm_alleles1):
        for j, alleles2 in enumerate(perm_alleles2):
            acc = sum(allele1.accuracy(allele2)[0] for (allele1, _), (allele2, _) in zip(alleles1, alleles2))
            if acc > best_acc:
                best_acc = acc
                best_i = i
                best_j = j
    return perm_alleles1[best_i], perm_alleles2[best_j]


def process(gene, sample, input_base, input_pred, avail, loo, out):
    if input_base and input_pred:
        base_alleles, pred_alleles = permute_best(input_base, input_pred)
    else:
        base_alleles = input_base or ()
        pred_alleles = input_pred or ()

    for base, pred in itertools.zip_longest(base_alleles, pred_alleles):
        if base is None:
            assert pred is not None
            base_allele = '<INPUT_MISSING>' if input_base is None else '<MISSING>'
            acc = 0
            base_len = np.nan
            curr_avail = 'NA' if avail is None else 'T'
        else:
            base_allele = base[0]
            if avail is None:
                curr_avail = 'NA'
            elif all(opt[0] is None for opt in base_allele.options):
                curr_avail = 'T'
            else:
                curr_avail = any(
                    not loo or sample != sample2
                        for opt in base_allele.options
                        for sample2 in avail.get(opt[:AVAIL_SIZE], ())
                    )
                curr_avail = 'T' if curr_avail else 'F'

        if pred is None:
            pred_allele = '<INPUT_MISSING>' if input_pred is None else '<MISSING>'
            qual = 0
            acc = 0
            base_len = np.nan
        else:
            pred_allele, qual = pred

        if base is not None and pred is not None:
            acc, base_len = base_allele.accuracy(pred_allele)
        out.write(f'{gene}\t{sample}\t{base_allele}\t{pred_allele}\t{qual}\t{acc}\t{base_len}\t{curr_avail}\n')


AVAIL_SIZE = 2


def get_available_alleles(baseline, avail_samples):
    avail = {}
    for (gene, sample), alleles in baseline.items():
        if avail_samples is not None and sample not in avail_samples:
            continue
        for allele, _ in alleles:
            for opt in allele.options:
                subopt = opt[:AVAIL_SIZE]
                # Safer than using defaultdict later.
                if gene not in avail:
                    avail[gene] = {}
                if subopt not in avail[gene]:
                    avail[gene][subopt] = []
                avail[gene][subopt].append(sample)
    return avail


def load_avail(filename):
    if filename is None:
        return None
    with common.open(filename) as f:
        avail_dict = load_alleles(f)
        avail = {}
        for (gene, sample), alleles in avail_dict.items():
            for allele, _ in alleles:
                for opt in allele.options:
                    subopt = opt[:AVAIL_SIZE]
                    # Safer than using defaultdict later.
                    if gene not in avail:
                        avail[gene] = {}
                    if subopt not in avail[gene]:
                        avail[gene][subopt] = []
                    avail[gene][subopt].append(sample)
        return avail


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--base', metavar='FILE', required=True,
        help='Baseline alleles.')
    parser.add_argument('-p', '--pred', metavar='FILE', required=True,
        help='Predicted alleles.')
    parser.add_argument('--loo', action='store_true',
        help='Consider leave-one-out for the set of available alleles.')
    parser.add_argument('--avail', metavar='FILE', required=False,
        help='Only provided samples were available during genotyping.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    args = parser.parse_args()

    with common.open(args.base) as f:
        base = load_alleles(f)
    with common.open(args.pred) as f:
        pred = load_alleles(f)

    avail = None
    if args.avail:
        with common.open(args.avail) as f:
            avail_samples = set(map(str.strip, f))
        avail = get_available_alleles(base, avail_samples)
    elif args.loo:
        avail = get_available_alleles(base, None)

    with common.open(args.output, 'w') as out:
        out.write('gene\tsample\tbase\tpred\tqual\taccuracy\tbase_length\tavail\n')
        for key in sorted(set(base.keys()) | set(pred.keys())):
            gene, sample = key
            process(gene, sample, base.get(key), pred.get(key), avail and avail.get(gene), args.loo, out)


if __name__ == '__main__':
    main()
