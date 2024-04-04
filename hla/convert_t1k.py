#!/usr/bin/env python3

import os
import sys
import glob
import argparse
from collections import defaultdict

import common


def common_prefix(s1, s2):
    for i, (x, y) in enumerate(zip(s1, s2)):
        if x != y:
            return i
    return min(len(s1), len(s2))


def process_alleles(s):
    if s == '.':
        return ('<MISSING>',)
    return [x.split('*', 1)[1] for x in s.split(',')]


def specify_alleles(f, gene_alleles):
    for line in f:
        gene, allele = line.strip().split()[0].split('*')
        alleles1, _, alleles2, _ = gene_alleles[gene]

        for i, allele1 in enumerate(alleles1):
            if allele.startswith(allele1):
                alleles1[i] = allele
                allele = None
                break
        if allele is None:
            continue
        for j, allele2 in enumerate(alleles2):
            if allele.startswith(allele2):
                alleles2[j] = allele
                allele = None
                break
        # assert allele is None, f'{f.name}: "{gene}*{allele}"'


def process_sample(dir, out):
    dir_split = dir.split('/')
    sample = dir_split[dir_split.index('.') + 1]
    filenames = glob.glob(os.path.join(dir, '*_genotype.tsv'))
    assert len(filenames) == 1

    gene_alleles = {}
    with common.open(filenames[0]) as f:
        for line in f:
            line = line.strip().split('\t')
            gene = line[0]

            alleles1 = process_alleles(line[2])
            qual1 = float(line[4])
            alleles2 = process_alleles(line[5])
            qual2 = float(line[7])

            if len(line) > 8:
                extra_alleles = line[8].split(';')[0]
                for allele3 in extra_alleles.split(','):
                    allele3 = allele3.split('*', 1)[1]
                    if max(common_prefix(allele3, x) for x in alleles1) \
                            >= max(common_prefix(allele3, y) for y in alleles2):
                        alleles1.append(allele3)
                    else:
                        alleles2.append(allele3)
            assert gene not in gene_alleles
            gene_alleles[gene] = (alleles1, qual1, alleles2, qual2)

    filenames = glob.glob(os.path.join(dir, '*_allele.tsv'))
    assert len(filenames) == 1
    with common.open(filenames[0]) as f:
        specify_alleles(f, gene_alleles)

    for gene, (alleles1, qual1, alleles2, qual2) in gene_alleles.items():
        out.write(f'{gene}\t{sample}\t')
        out.write('{} {:.3g}\t{} {:.3g}\n'.format('|'.join(alleles1), qual1, '|'.join(alleles2), qual2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--t1k', required=True, nargs='+',
        help='Multiple T1K output folders, one for each sample '
            '(sample identified by `./` before it).')
    parser.add_argument('-o', '--output', required=True,
        help='Converted T1K predictions.')
    args = parser.parse_args()

    res = []
    with common.open(args.output, 'w') as out:
        out.write(f'# {" ".join(sys.argv)}\n')
        out.write('gene\tsample\talleles\n')
        for dir in args.t1k:
            res.append(process_sample(dir, out))


if __name__ == '__main__':
    main()
