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

## local:

import intypes

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
    "Assigns new gene_ids based on p_exon_overlap and "
    "p_exons_overlap. When a transcript is found in both input "
    "files, but the exon boundaries differ (within the specified "
    "tolerances), the output exon coordinates will be taken from "
    "primary_gtf. Also outputs a tab-delimited table "
    "cross-referencing filtered transcript_ids to the kept "
    "transcript_id with which they were found to be redundant."
)

def build_parser():

    parser = argparse.ArgumentParser(
	description = description,
	epilog = epilog,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
	"primary_gtf", 
	help="GTF2.2 file containing primary isoform list"
    )

    parser.add_argument(
	"secondary_gtf", 
	nargs='?', 
	help="GTF2.2 file with secondary isoform list"
    )

    parser.add_argument(
	"--tol_sj", 
	type=intypes.non_negative_int, 
	default=0,
	help="Tolerance (bp) for matching splice junction coordinates"
    )

    parser.add_argument(
	"--tol_tss", 
	type=intypes.non_negative_int, 
	default=0,
	help="Tolerance (bp) for matching transcript start coordinates"
    )

    parser.add_argument(
	"--tol_tts", 
	type=intypes.non_negative_int, 
	default=0,
	help="Tolerance (bp) for matching transcript end coordinates"
    )

    parser.add_argument(
	"--p_exon_overlap",
	type=intypes.float_proportion,
	default=0.5,
	help="Minimum proportion overlap between two exons for gene matching"
    )

    parser.add_argument(
	"--p_exons_overlap",
	type=intypes.float_proportion,
	default=0.2,
	help="Minimum proportion of exons with overlaps needed for gene matching"
    )

    parser.add_argument(
	"--keep_all_primary",
	action='store_true',
	help="Keep all primary transcripts (even redundant ones)"
    )

    parser.add_argument(
	"--gene_prefix",
	type=intypes.non_whitespace_str,
	default='LOC.',
	help="Prefix for gene_ids and transcript_ids"
    )

    parser.add_argument(
	"--primary_prefix",
	type=intypes.non_whitespace_str,
	default=None,
	help="Prefix to be prepended to primary_gtf transcript_ids"
    )

    parser.add_argument(
       "--secondary_prefix",
       type=intypes.non_whitespace_str,
       default=None,
       help="Prefix to be prepended to secondary_gtf transcript_ids"
    )

    parser.add_argument(
	"--output_prefix",
	type=intypes.non_whitespace_str,
	default='union',
	help="Prefix for output files"
    )

    parser.add_argument(
	"--ids_from",
	type=int,
	choices=[1, 2],
	default=1,
	help=(
	    "Source of transcript_ids and gene_ids for redundant transcripts; "
	    "1: primary_gtf; 2: secondary_gtf"
	)
    )

    return parser


def parse_args(parser): 

    args = parser.parse_args()

    arg_names = [
        'primary_gtf',
        'secondary_gtf',
        'tol_sj',
        'tol_tss',
        'tol_tts',
        'keep_all_primary',
        'p_exon_overlap',
        'p_exons_overlap',
        'primary_prefix',
        'secondary_prefix',
        'output_prefix',
        'ids_from', 
    ]

    for arg_name in arg_names:
        print(f"{arg_name}: {getattr(args, arg_name)}")

    return args


def initialize():

    time_start = time.time()

    parser = build_parser()
    
    params = parse_args(parser)
    params.time_start = time_start

    time_stamp = time.strftime(
        '%Y%m%d%H%M%S',
        time.localtime(time_start)
    )

    print(
        f"time_stamp: {time_stamp}\n"
        f"hostname: {socket.gethostname()}\n"
        f"working_directory: {os.getcwd()}"
    )

    return params


if __name__ == "__main__":
    params = initialize()
    print(f"initialize() yielded: {params}")
    sys.exit(0)

