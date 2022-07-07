import argparse
import os
import socket
import sys
import time

description = (
    "Merges primary_gtf with secondary_gtf, resulting in a non-redundant union gtf. "
    "Omitting secondary_gtf results in non-redundant set from primary_gtf."
)

epilog = (
    "Output includes 'exon' and 'transcript' features, "
    "with attributes (in order) "
    "'gene_id', 'transcript_id', 'primary_id', 'secondary_id; "
    "primary_id is the transcript_id in primary_gtf, and "
    "secondary_id is the transcript_id in secondary_gtf; "
    "set to empty string '' if not present in corresponding gtf. "
    "Error if transcript_id in primary_gtf is also "
    "found in secondary_gtf; set --primary_prefix and/or "
    "--secondary_prefix to fix this."
)

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
    type=int, 
    default=0,
    help="Tolerance (bp) for matching splice junction coordinates"
)

parser.add_argument(
    "--tol_tss", 
    type=int, 
    default=0,
    help="Tolerance (bp) for matching transcriptional start coordinates"
)

parser.add_argument(
    "--tol_tts", 
    type=int, 
    default=0,
    help="Tolerance (bp) for matching transcript end coordinates"
)

parser.add_argument(
    "--ids_from", 
    type=int, 
    choices=[1, 2], 
    default=1,
    help="Where to get transcript_ids and gene_ids; 1: primary_gtf; 2: secondary_gtf"
)

parser.add_argument(
    "--keep_all_primary",
    action='store_true',
    help="Keep all primary transcripts (even redundant ones)"
)

parser.add_argument(
    "--primary_prefix",
    type=str,
    default=None,
    help="Prefix for primary_gtf transcript_ids"
)

parser.add_argument(
   "--secondary_prefix",
   type=str,
   default=None,
   help="Prefix for secondary_gtf transcript_ids"
)

def initialize():

    time_start = time.time()

    time_stamp = time.strftime(
        '%Y%m%d%H%M%S',
        time.localtime(time_start)
    )

    sys.stderr.write(
        f"time_stamp: {time_stamp}\n"
        f"hostname: {socket.gethostname()}\n"
        f"working_directory: {os.getcwd()}\n"
    )

    return time_start


def parse_args(): 

    args = parser.parse_args()

    arg_names = [
        'primary_gtf',
        'secondary_gtf',
        'tol_sj',
        'tol_tss',
        'tol_tts',
        'ids_from',
        'keep_all_primary',
        'primary_prefix',
        'secondary_prefix'
    ]

    for arg_name in arg_names:
        sys.stderr.write(f"{arg_name}: {getattr(args, arg_name)}\n")

    return args


if __name__ == "__main__":
    initialize()
    args = parse_args()
    print(f"parse_args() yielded: {args}")
    sys.exit(0)

