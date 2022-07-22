#!/usr/bin/env python
# cython: embedsignature=True, binding=True

"""
Contains initialization and argument handling code. 
"""

import argparse
import os
import re
import socket
import sys
import time

def float_proportion(val):

    try:
        val = float(val)
        if not (0 <= val <= 1):
            raise Exception()
    except Exception as e:
        raise argparse.ArgumentTypeError(
            f"val '{val}' not a proportion in the interval [0, 1]"
        )

    return val


def strictly_positive_int(val):

    try:
        val = int(val)
        if val < 1:
            raise Exception()
    except Exception as e:
        raise argparse.ArgumentTypeError(
            f"val '{val}' not a strictly positive integer"
        )    

    return val


def non_negative_int(val):

    try:
        val = int(val)
        if val < 0:
            raise Exception()
    except Exception as e:
        raise argparse.ArgumentTypeError(
            f"val '{val}' not a non-negative integer"
        )

    return val


def non_whitespace_str(val):

    pat = re.compile(r'\s')

    try:
        val = str(val)
        if pat.search(val):
            raise Exception()
    except Exception as e:
        raise argparse.ArgumentTypeError(
            f"val '{val}' should be a string with no whitespace"
        )

    return val
    

description = (
    "Merges redundant transcripts from GTF2.2 formatted files. "
    "Merges primary_gtf with secondary_gtf, resulting in a non-redundant union gtf. "
    "Omitting secondary_gtf results in non-redundant set from primary_gtf."
)

epilog = (
    "Input GTF2.2 formatted files are required to have exon "
    "features with attributes that include 'gene_id' and "
    "'transcript_id'. Output GTF2.2 formatted file includes "
    "'exon' and 'transcript' features, both with attributes "
    "(in order) 'gene_id', 'transcript_id', 'primary_id', "
    "'secondary_id; where primary_id is the transcript_id "
    "in primary_gtf, and secondary_id is the transcript_id "
    "in secondary_gtf; primary_id and/or secondary_id will be "
    "set to empty string '' if the transcript is not present "
    "in the corresponding GTF file. Throws error if transcript_id in "
    "primary_gtf is also found in secondary_gtf; set "
    "--primary_prefix and/or --secondary_prefix to fix this. "
    "When a transcript is found in both input files, but the "
    "exon boundaries differ (within the specified tolerances), "
    "the output exon coordinates will be taken from "
    "primary_gtf. Also outputs a tab-delimited table "
    "cross-referencing filtered transcript_ids to the kept "
    "transcript_id with which they were found to be redundant."
)

parser = argparse.ArgumentParser(
    description = description,
    epilog = epilog,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "gtf_list_file", 
    help="File with list of GTF2.2 files (one path per line) to be merged"
)

parser.add_argument(
    '--nthreads',
    type=strictly_positive_int,
    default=1,
    help="Number of threads to use for parallel execution; even number more efficient"
)

parser.add_argument(
    '--memory_gb',
    type=strictly_positive_int,
    default=8,
    help="Amount of available RAM, in gigabytes"
)

parser.add_argument(
    "--tol_sj", 
    type=non_negative_int, 
    default=0,
    help="Tolerance (bp) for matching splice junction coordinates"
)

parser.add_argument(
    "--tol_tss", 
    type=non_negative_int, 
    default=0,
    help="Tolerance (bp) for matching transcript start coordinates"
)

parser.add_argument(
    "--tol_tts", 
    type=non_negative_int, 
    default=0,
    help="Tolerance (bp) for matching transcript end coordinates"
)

parser.add_argument(
    "--p_exon_overlap",
    type=float_proportion,
    default=0.5,
    help="Minimum proportion overlap between two exons for gene matching"
)

parser.add_argument(
    "--p_exons_overlap",
    type=float_proportion,
    default=0.25,
    help="Minimum proportion of exons with overlaps needed for gene matching"
)

parser.add_argument(
    "--output_prefix",
    type=non_whitespace_str,
    default='union',
    help="Prefix for output files"
)


def parse_args(): 

    args = parser.parse_args()

    arg_names = [
        'gtf_list_file',
        'nthreads',
        'memory_gb',
        'tol_sj',
        'tol_tss',
        'tol_tts',
        'p_exon_overlap',
        'p_exons_overlap',
        'output_prefix',
    ]

    for arg_name in arg_names:
        val = getattr(args, arg_name)
        print(f"{arg_name}: {val}")

    return args


def initialize():

    time_start = time.time()

    time_stamp = time.strftime(
        '%Y%m%d%H%M%S',
        time.localtime(time_start)
    )

    params = parse_args()
    params.time_start = time_start
    params.memory = params.memory_gb * 10**9

    print(
        f"time_stamp: {time_stamp}\n"
        f"hostname: {socket.gethostname()}\n"
        f"working_directory: {os.getcwd()}\n"
        f"memory: {params.memory}\n"
    )

    return params


if __name__ == "__main__":
    initialize()
    args = parse_args()
    print(f"parse_args() yielded: {args}")
    sys.exit(0)

