#!/usr/bin/env python3

import re
import argparse
import collections
import os

import common


def parse_annotation(f):
    res = {}
    for line in f:
        if line.startswith('## gene (copy num = 0):'):
            for gene in line.strip().split(':')[1].split(','):
                res[gene.strip()] = ('<MISSING>',)
        elif line.startswith('## gene (copy num'):
            for gene in line.strip().split(':')[1].split(','):
                res[gene.strip()] = []

        if line.startswith('#'):
            continue

        line = line.strip().split('\t')
        if line[2] != 'transcript':
            continue
        annot = dict(re.findall(r'([a-zA-Z0-9_-]+)\s+"([^"]+)"\s*;', line[8]))
        if 'consensus' not in annot:
            continue

        gene = annot['gene_name']
        assert ',' not in annot['consensus']
        alleles = '|'.join(allele.strip().split('*', 1)[1] for allele in annot['alleles'].split(','))
        if gene not in res:
            res[gene] = []
        res[gene].append(alleles)
    return res


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--annotation', metavar='FILE', required=True, nargs='+',
        help='Multiple annotation files with prefixes {sample}.{hap}.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output filename.')
    args = parser.parse_args()

    annotations = collections.defaultdict(lambda: [None, None])
    for filename in args.annotation:
        name_split = os.path.basename(filename).split('.')
        sample = name_split[0]
        try:
            hap = int(name_split[1])
        except ValueError:
            hap = 0
        print(filename)
        with common.open(filename) as f:
            annotations[sample][max(0, hap - 1)] = parse_annotation(f)

    with open(args.output, 'w') as out:
        out.write('gene\tsample\talleles1\talleles2\n')
        for sample, (annot1, annot2) in sorted(annotations.items()):
            for gene in sorted(set(annot1.keys()) | (set(annot2.keys()) if annot2 else set())):
                alleles1 = ';'.join(annot1.get(gene) or ('<UNKNOWN>',))
                alleles2 = ';'.join(annot2.get(gene) or ('<UNKNOWN>',)) if annot2 else '<HAPLOID>'
                out.write(f'{gene}\t{sample}\t{alleles1}\t{alleles2}\n')


if __name__ == '__main__':
    main()
