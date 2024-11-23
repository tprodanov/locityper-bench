#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import random

import common


OPTIONS = {
    'A': 'CGT',
    'C': 'AGT',
    'G': 'ACT',
    'T': 'ACG',
}


def mutate(seq, rate, always):
    s = ''
    count = 0
    for nt in seq:
        if nt != 'N' and random.random() < rate:
            s += random.choice(OPTIONS[nt])
            count += 1
        else:
            s += nt
    if count or not always:
        return s, count

    while True:
        i = random.randrange(len(seq))
        nt = seq[i]
        if nt != 'N':
            new = random.choice(OPTIONS[nt])
            return s[:i] + new + s[i+1:], 1


def write_fa(name, seq, out):
    out.write(f'>{name}\n')
    for i in range(0, len(seq), 120):
        out.write(seq[i : i + 120] + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='Input FASTA file.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output FASTA file.')
    parser.add_argument('-m', '--mut-rate', metavar='FLOAT', default=0.0001, type=float,
        help='Chance for every given nucleotide to mutate in a given haplotype.')
    parser.add_argument('-c', '--copies', metavar='INT', default=1, type=int,
        help='Total number of mutated copies (in total) per haplotype.')
    parser.add_argument('-s', '--seed', metavar='INT', type=int,
        help='Optional: random seed.')
    parser.add_argument('-a', '--always', action='store_true',
        help='If no mutations were added, add one randomly.')
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    with common.open(args.input) as inp, common.open(args.output, 'w') as out:
        for record in SeqIO.parse(inp, 'fasta'):
            name = record.id
            seq = str(record.seq)
            for i in range(1, args.copies + 1):
                new_seq, count = mutate(seq, args.mut_rate, args.always)
                write_fa(f'{name}-mut{i} {count} mutations', new_seq, out)


if __name__ == '__main__':
    main()
