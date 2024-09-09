#!/usr/bin/env python3

import os
import sys
import multiprocessing
import time
import numpy as np
from collections import defaultdict
import pysam
import argparse


def find_bams(d):
    for f in os.listdir(d):
        subdir = os.path.join(d, f)
        bam_filename = os.path.join(subdir, 'alns', '00.bam')
        success_filename = os.path.join(subdir, 'success')
        if os.path.exists(bam_filename) and os.path.exists(success_filename):
            yield (f, bam_filename)
        else:
            sys.stderr.write(f'ERROR: Incomplete directory {subdir}\n')


def find_all_bams(d, samples):
    filenames = []
    for sample in map(str.strip, samples):
        for locus, filename in find_bams(os.path.join(d, sample, 'loci')):
            filenames.append((sample, locus, filename))
    return filenames


def process_one(sample, locus, filename):
    try:
        with pysam.AlignmentFile(filename) as bam:
            contigs = { contig: i for i, contig in enumerate(sorted(bam.header.references)) }
            contigs[None] = 2
            locs = defaultdict(lambda: [0, 0, 0])
            locs2 = defaultdict(lambda: [-np.inf, -np.inf])
            for record in bam:
                loc = contigs[record.reference_name]
                if record.has_tag('pr'):
                    pr = float(record.get_tag('pr'))
                    locs[record.query_name][loc] += pr
                else:
                    aln_lik = float(record.get_tag('al'))
                    locs2[record.query_name][loc] = max(locs2[record.query_name][loc], aln_lik)

        s = ''
        for name, (pr1, pr2, pr3) in locs.items():
            assert name not in locs2
            s += f'{sample}\t{locus}\t{name}\t{pr1:0.2f}\t{pr2:0.2f}\t{pr3:0.2f}\n'
        for name, (lik1, lik2) in locs2.items():
            diff = lik1 - lik2
            if diff > 1:
                pr1, pr2 = (0.8, 0.2)
            elif diff < -1:
                pr1, pr2 = (0.2, 0.8)
            else:
                pr1, pr2 = (0.5, 0.5)
            s += f'{sample}\t{locus}\t{name}\t{pr1:0.2f}\t{pr2:0.2f}\tNA\n'

        return s
    except Exception as e:
        sys.stderr.write(f'Encountered error while processing {sample} at {locus} ({filename}):\n{e}\n')
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='DIR',
        help='Input directory, which contains analysis for multiple samples.')
    parser.add_argument('-s', '--samples', metavar='FILE', required=True,
        help='File with samples to process.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output file.')
    parser.add_argument('-@', '--threads', metavar='INT', type=int, default=8,
        help='Number of processing threads [%(default)s].')
    args = parser.parse_args()

    with open(args.samples) as f:
        all_paths = find_all_bams(args.input, f)

    out = open(args.output, 'w')
    out.write('sample\tlocus\tread_name\th1\th2\tunmapped\n')

    DELTA = 10.0
    time_thresh = time.perf_counter() + DELTA
    finished = 0
    errors = 0
    total = len(all_paths)
    total_len = len(f'{total:,}')

    def callback(s):
        nonlocal finished, errors, time_thresh
        finished += 1
        if s is None:
            errors += 1
        else:
            out.write(s)

        curr_time = time.perf_counter()
        if time_thresh < curr_time:
            time_thresh = curr_time + DELTA
            sys.stderr.write(f'[{finished:{total_len},} / {total:,}]\n')
            out.flush()

    with multiprocessing.Pool(args.threads) as pool:
        results = [pool.apply_async(process_one, tup, callback=callback) for tup in all_paths]
        for res in results:
            res.get()
        pool.close()
        pool.join()
    sys.stderr.write(f'Successfully processed {finished - errors} sample-loci pairs\n')
    if errors:
        sys.stderr.write(f'Encountered errors in {errors} pairs\n')


if __name__ == '__main__':
    main()
