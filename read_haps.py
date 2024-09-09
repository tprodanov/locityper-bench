#!/usr/bin/env python3

import os
import sys
import multiprocessing
import time
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
    with pysam.AlignmentFile(filename) as bam:
        contigs = { contig: i for i, contig in enumerate(sorted(bam.header.references)) }
        contigs[None] = 2
        locs = defaultdict(lambda: [0, 0, 0])
        for record in bam:
            loc = contigs[record.reference_name]
            pr = float(record.get_tag('pr'))
            locs[record.query_name][loc] += pr

    s = ''
    for name, (pr0, pr1, pr2) in locs.items():
        s += f'{sample}\t{locus}\t{name}\t{pr0:0.2f}\t{pr1:0.2f}\t{pr2:0.2f}\n'
    return s


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
    total = len(all_paths)

    def callback(s):
        nonlocal finished, time_thresh
        finished += 1
        out.write(s)

        curr_time = time.perf_counter()
        if time_thresh < curr_time:
            time_thresh = curr_time + DELTA
            sys.stderr.write(f'Finished [{finished:4} / {total}]\n')
            out.flush()

    with multiprocessing.Pool(args.threads) as pool:
        results = [pool.apply_async(process_one, tup, callback=callback) for tup in all_paths]
        for res in results:
            res.get()
        pool.close()
        pool.join()
    sys.stderr.write(f'Successfully processed {total} sample-loci pairs\n')


if __name__ == '__main__':
    main()
