#!/usr/bin/env python

"""Convert PQC text format to JSON format.

Synopsis

  python txt2json.py input.txt > output.json

"""

import argparse
import json
import os
import re
import sys

ROW_SEPERATOR = '\t'

META_REGEX = re.compile(r'([\w]+)\s*\:\s*(.*)')
HEADER_REGEX = re.compile('([\w]+)(?:\[([^\]]+)\])?')
BOOL_PATTERN = {'false': False, 'true': True}

def split_row(line):
    return line.split(ROW_SEPERATOR)

def parse_list(value, type=str, separator=','):
    """Parse list like expression.
    >>> parse_list("(2, 4, 7)", type=float)
    [2.0, 4.0, 7.0]
    """
    return [type(pos.strip()) for pos in value.strip().strip('()[]').split(separator)]

def fetch_meta(meta, series_units, series, line):
    result = META_REGEX.match(line)
    if result:
        key, value = result.groups()
        if key in meta:
            raise KeyError(f"key already exists: '{key}'")
        try:
            value = float(value)
            if value.is_integer():
                value = int(value)
        except ValueError:
            value = BOOL_PATTERN.get(value.lower(), value)
        if key == 'table_position':
            value = parse_list(value, type=float)
        meta[key] = value
        return fetch_meta
    return fetch_header(meta, series_units, series, line)

def fetch_header(meta, series_units, series, line):
    for item in split_row(line):
        item = item.strip()
        result = HEADER_REGEX.match(item)
        if not result:
            raise ValueError(f"invalid header: {items}")
        key, value = result.groups()
        series_units[key] = value
        series[key] = []
    return fetch_row

def fetch_row(meta, series_units, series, line):
    items = split_row(line)
    row = []
    for i, key in enumerate(series.keys()):
        value = float(items[i])
        series.get(key).append(value)
        row.append(value)
    return fetch_row

def load_text(f):
    meta = {}
    series_units = {}
    series = {}
    func = fetch_meta
    for line in f:
        line = line.strip()
        func = func(meta, series_units, series, line)
    return {'meta': meta, 'series_units': series_units, 'series': series}

def to_json(data, f):
    json.dump(data, f, indent=2)
    sys.stdout.write(os.linesep)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=argparse.FileType('r'), help="Input file in PQC text format.")
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help="Write output to file (default is stdout).")
    return parser.parse_args()

def main():
    args = parse_args()
    with args.file as f:
        data = load_text(f)
    with args.o as f:
        to_json(data, f)

if __name__ == '__main__':
    main()
