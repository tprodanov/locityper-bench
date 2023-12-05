#!/usr/bin/env python3

import argparse
import pysam
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import common


def add_fasta_source(source, filename, out_seqs):
    with common.open(filename) as fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            out_seqs.append(SeqRecord(
                record.seq,
                id=f'{source}::{record.id}',
                description=f'{len(record.seq)} bp'
            ))


def add_vcf_source(source, filename, ref_seq, chrom, start, end, out_seqs):
    # Output indels and SNPs over 5 bp.
    SNP_LEN = 5

    with pysam.VariantFile(filename) as vcf:
        prev_end = 0
        for var in vcf.fetch(chrom, start, end):
            assert var.start >= prev_end, f'Source {source} ({filename}) contains overlapping variants'
            alleles = var.alleles
            ref_allele = alleles[0]
            if max(map(len, alleles)) < SNP_LEN and all(len(allele) == len(ref_allele) for allele in alleles[1:]):
                continue

            var_end = var.start + len(ref_allele)
            assert start < var_end and var.start < end
            prefix = ref_seq[ : max(0, var.start - start)]
            suffix = ref_seq[min(len(ref_seq), var_end - start) : ]
            for i, alt in enumerate(alleles[1:], 1):
                seq = prefix + alt + suffix
                out_seqs.append(SeqRecord(
                    Seq(seq),
                    id=f'{source}::{var.start + 1}:{i}',
                    description=f'{len(seq)} bp, variant {chrom}:{var.start + 1}-{var_end} ({len(ref_allele)} bp), '
                        f'{len(alleles)} alleles'))


def main():
    parser = argparse.ArgumentParser(description='Extract relevant sequences from various sources.')
    parser.add_argument('-r', '--reference', metavar=('NAME', 'FILE'), nargs=2, required=True,
        help='Reference name and FASTA file.')
    parser.add_argument('-l', '--locus', metavar='FILE', required=True,
        help='Locus coordinates in BED format.')
    parser.add_argument('-s', '--source', metavar='FILE', action='append',
        help='Source in FASTA or VCF.GZ format. '
            'Fasta sequence should be already reverse-complemented, if needed.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output FASTA file.')
    args = parser.parse_args()

    with common.open(args.locus) as f:
        line = next(f).strip().split()
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])

    ref = pysam.FastaFile(args.reference[1])
    ref_seq = ref.fetch(chrom, start, end).upper()
    out_seqs = [SeqRecord(
        Seq(ref_seq),
        id=args.reference[0],
        description=f'{chrom}:{start+1}-{end}, {len(ref_seq)} bp')]

    sources = set()
    for filename in args.source:
        source = os.path.basename(filename).split('.', 1)[0]
        if source in sources:
            assert False, f'Source {source} appears twice'
        sources.add(source)

        if filename.endswith('.fa') or filename.endswith('.fasta') \
                or filename.endswith('.fa.gz') or filename.endswith('.fasta.gz'):
            add_fasta_source(source, filename, out_seqs)
        else:
            add_vcf_source(source, filename, ref_seq, chrom, start, end, out_seqs)

    with common.open(args.output, 'w') as out:
        SeqIO.write(out_seqs, out, 'fasta')


if __name__ == '__main__':
    main()
