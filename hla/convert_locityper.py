#!/usr/bin/env python3

import argparse
import collections
import sys

import common


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='Locityper CSV file for multiple loci and samples.')
    parser.add_argument('-a', '--annotation', metavar='FILE', required=True,
        help='CSV file with annotation for database samples.')
    parser.add_argument('-l', '--loci', metavar='FILE', required=True,
        help='Two column file with gene-locus correspondence.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    args = parser.parse_args()

    with common.open(args.loci) as f:
        genes = []
        locus_to_genes = collections.defaultdict(list)
        for line in f:
            gene, locus = line.strip().split('\t')
            genes.append(gene)
            locus_to_genes[locus].append(gene)

    with common.open(args.annotation) as f:
        true_gts = { gene: {} for gene in genes }
        for row in common.read_csv(f):
            curr_gts = true_gts.get(row['gene'])
            if curr_gts is not None:
                curr_gts[row['sample']] = (row['alleles1'].split(';'), row['alleles2'].split(';'))

    with common.open(args.input) as inp, common.open(args.output, 'w') as out:
        out.write(f'# {" ".join(sys.argv)}\n')
        out.write('gene\tsample\talleles\n')
        for row in common.read_csv(inp):
            locus = row['locus']
            if locus not in locus_to_genes:
                continue

            sample = row['sample']
            gt = row['genotype']
            warnings = row['warnings']
            weight_dist = float(row['weight_dist'])
            total_reads = float(row['total_reads'])
            unexpl_reads = float(row['unexpl_reads'])
            if warnings == '*' and weight_dist <= 20 and (unexpl_reads < 1000 or unexpl_reads <= 0.1 * total_reads):
                qual = 1
            else:
                qual = 0

            for gene in locus_to_genes[locus]:
                curr_alleles = []
                gene_true_gts = true_gts[gene]
                for allele in gt.split(','):
                    if allele in gene_true_gts:
                        curr_alleles.extend(gene_true_gts[allele][0])
                    elif '.' in allele:
                        pred_sample, hap = allele.split('.')
                        curr_alleles.extend(gene_true_gts[pred_sample][int(hap) - 1])
                    else:
                        if '*' in allele:
                            allele = allele.split('*', 1)[1]
                        curr_alleles.append(allele)

                out.write(f'{gene}\t{sample}')
                for allele in curr_alleles:
                    out.write(f'\t{allele} {qual}')
                out.write('\n')


if __name__ == '__main__':
    main()
