import sys
import gzip
import os
import builtins
import csv


def open(filename, mode='r'):
    assert mode == 'r' or mode == 'w'
    if filename is None or filename == '-':
        return sys.stdin if mode == 'r' else sys.stdout
    elif filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    else:
        return builtins.open(filename, mode)


def read_csv(f):
    """
    Wrapper over csv.DictReader, that ignores first lines with comments.
    """
    fields = None
    for line in f:
        if not line.startswith('#'):
            fields = line.strip().split('\t')
            break
    assert fields is not None
    return csv.DictReader(f, fieldnames=fields, delimiter='\t')

