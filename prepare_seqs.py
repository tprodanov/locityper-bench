#!/usr/bin/env python3

import os
import sys
import operator
from collections import OrderedDict, defaultdict
import argparse
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import Seq

import common


def get_params(record, alns, qseqs):
    p = argparse.Namespace()
    p.qname = record.query_name
    p.qstart = record.query_alignment_start
    p.qend = record.query_alignment_end

    if '.' in p.qname:
        p.sample, hap = p.qname.split('.')
        p.hap = int(hap) - 1
    else:
        p.sample = p.qname
        p.hap = 0

    p.rname = record.reference_name
    p.rstart = record.reference_start
    p.rend = record.reference_end
    p.rlen = alns.get_reference_length(p.rname)
    p.edit = record.get_tag('NM')
    p.nmatches = p.rend - p.rstart - p.edit
    p.rrem = p.rlen - p.rend
    p.clip = (p.rstart + p.rrem) / p.rlen
    p.divergence = (p.edit + p.rstart + p.rrem) / p.rlen

    if p.qname not in qseqs:
        assert record.seq is not None and 'H' not in record.cigarstring
        qseqs[p.qname] = record.seq
    p.qseq = qseqs[p.qname]
    p.qlen = len(p.qseq)
    return p


def process_alns(alns, seqs, out, clipping, max_divergence, add_padding):
    coords = defaultdict(list)
    gts = defaultdict(lambda: (set(), set()))
    locus = next(iter(seqs.values())).split('*')[0]
    n_ext = 0
    qseqs = {}
    for record in alns:
        p = get_params(record, alns, qseqs)
        if record.is_unmapped:
            gts[p.sample][p.hap].add('<MISSING>')
            sys.stderr.write(f'!!! {p.qname} is unmapped\n')
            out.write(f'{p.sample}\t{p.hap+1}\t*\t*\t*\tNA\tNA\tNA\tNA\t!!\t*\n')
            continue

        if p.qname not in coords and p.qname == p.sample:
            gts[p.sample][1].add('<HAPLOID>')
        is_redundant = False
        for prev_start, prev_end, prev_matches in coords[p.qname]:
            overl_frac = max(0, min(prev_end, p.qend) - max(prev_start, p.qstart)) / p.rlen
            # Alignment is not good enough or overlaps an old one.
            if p.nmatches + 10 <= prev_matches or overl_frac >= 0.5:
                is_redundant = True
        if is_redundant:
            continue

        coords[p.qname].append((p.qstart, p.qend, p.nmatches))
        allele = None
        status = ''
        if p.clip >= clipping:
            status += '!'
            allele = '<MISSING>'
        elif p.clip >= 0.005:
            status += '~'
        else:
            status += '+'

        if p.divergence >= max_divergence:
            status += '!'
            allele = '<MISSING>'
        elif p.divergence >= 0.01:
            status += '~'
        else:
            status += '+'

        if allele is None:
            if record.seq is None:
                assert 'H' not in record.cigarstring
                assert not record.is_reverse
                seq = p.qseq
            else:
                seq = record.seq

            if add_padding:
                assert not record.is_reverse
                seq = seq[max(0, p.qstart - p.rstart):min(p.qlen, p.qend + p.rrem)]
            else:
                assert p.qend <= len(seq)
                seq = seq[p.qstart:p.qend]
            assert seq

            allele = seqs.get(seq)
            if allele is None:
                n_ext += 1
                allele = f'{locus}*ext{n_ext:03}'
                assert allele not in seqs
                seqs[seq] = allele
        gts[p.sample][p.hap].add(allele)
        out.write(f'{p.sample}\t{p.hap+1}\t{p.qstart+1}-{p.qend}\t{"-" if record.is_reverse else "+"}\t'
            f'{p.rname}\t{p.rstart},{p.rrem}\t{p.clip:.5g}\t'
            f'{p.edit}\t{p.divergence:.5g}\t{status}\t{allele}\n')
    return gts


def alleles_to_str(alleles):
    if not alleles:
        return '<UNKNOWN>'
    count_special = sum(allele.startswith('<') for allele in alleles)
    if count_special == len(alleles) or count_special == 0:
        return ';'.join(alleles)
    return ';'.join(allele for allele in alleles if not allele.startswith('<'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples', required=True,
        help='Full set of samples.')
    parser.add_argument('-r', '--ref', required=True,
        help='Reference sequences in FASTA format.')
    parser.add_argument('-a', '--aln', required=True,
        help='Input alignment file between assembly and IPD sequences.')
    parser.add_argument('-o', '--output', required=True,
        help='Output directory.')
    parser.add_argument('-c', '--clipping', type=float, default=0.01,
        help='Maximum allowed clipping as fraction of allele length [%(default)s].')
    parser.add_argument('-d', '--divergence', type=float, default=0.05,
        help='Maximum allowed divergence [%(default)s].')
    parser.add_argument('--add-padding', action='store_true',
        help='Add padding to new sequences, according to unmapped first/last several bp.')
    args = parser.parse_args()

    with common.open(args.ref) as inp:
        seqs = OrderedDict((str(record.seq).upper(), record.id) for record in SeqIO.parse(inp, 'fasta'))
    pysam_verbosity = pysam.set_verbosity(0)
    with pysam.AlignmentFile(args.aln, require_index=False) as alns, \
            common.open(os.path.join(args.output, 'alns.csv'), 'w') as out_alns:
        pysam.set_verbosity(pysam_verbosity)
        out_alns.write('sample\thaplotype\tcoords\tstrand\tallele'
            '\tclip\tfrac_clip\tedit\tdivergence\tstatus\tnew_allele\n')
        gts = process_alns(alns, seqs, out_alns, args.clipping, args.divergence, args.add_padding)

    with common.open(args.samples) as samples_f, \
            common.open(os.path.join(args.output, 'extended.fa.gz'), 'w') as out_fa, \
            common.open(os.path.join(args.output, 'true_gts.csv'), 'w') as out_gt:
        out_gt.write('sample\talleles1\talleles2\n')
        samples = list(map(str.strip, samples_f))
        for sample in sorted(set(samples) | set(gts)):
            alleles1 = alleles_to_str(gts[sample][0])
            alleles2 = alleles_to_str(gts[sample][1])
            out_gt.write(f'{sample}\t{alleles1}\t{alleles2}\n')
        for seq, name in seqs.items():
            out_fa.write(f'>{name} {len(seq)} bp\n{seq}\n')


if __name__ == '__main__':
    main()
