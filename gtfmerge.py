#!/usr/bin/env python3

## system:
import sys
import time

## local:
import initial
import gtf
import output
import resolve

"""
Description: Merges primary_gtf with secondary_gtf, resulting in a 
  non-redundant union gtf. Omitting secondary_gtf results in 
  non-redundant set from primary_gtf.

For usage: python3 gtfmerge.py -h

Internals: 
  0) gtfs must have exon features w/ annotation including gene_id and transcript_id
  1) get rid of duplicates in primary_gtf; record { filtered_id: kept_id, ... }
  2) get rid of duplicates in secondary_gtf; record { filtered_id: kept_id, ...}
  3) add all unique from primary_gtf to output
  4) find overlap between unique from primary_gtf and unique from secondary_gtf
  5) add all unique non-overlapping from secondary_gtf to output
  6) make table of cross-references/mappings and print to file
  7) print to file output gtf, assigning gene_id/transcript_id per id_from argument

  proper output order assumes locs and transcripts[transcript_id][5] [exon_keys] 
    stay sorted by (start, end);

Author: mitch.kostich@jax.org
"""

args = initial.parse_args()
time_start = initial.initialize()


################################################
## parse primary_gtf:

sys.stderr.write("\n" f"Reading primary_gtf '{args.primary_gtf}':\n")

try:
    with open(args.primary_gtf, 'r') as fh:
        primary_dat = gtf.parse_gtf(fh, args.primary_prefix)
except Exception as e:
    sys.stderr.write(
        f"ERROR: while processing primary_gtf: '{args.primary_gtf}':\n" 
        f"{e}\n"
    )
    sys.exit(31)

sys.stderr.write(
    f"  done at {round(time.time() - time_start, 2)} seconds\n"
    f"  {len(primary_dat['transcripts'].keys())} transcript_ids\n"
    f"  {len(primary_dat['exons'].keys())} exons\n"
    f"  {len(primary_dat['locs'].keys())} seqids\n"
    f"  max_loc_length: {primary_dat['max_loc_length']}\n\n"
)


################################################
## resolve primary_gtf against itself:

if not args.keep_all_primary:
    sys.stderr.write("Resolving primary_gtf:\n")
    primary_dat = resolve.resolve_self(primary_dat, 1, args)
    sys.stderr.write(
        f"  done at {round(time.time() - time_start, 2)} seconds\n"
        f"  {len(primary_dat['transcripts'].keys())} transcript_ids\n"
        f"  {len(primary_dat['exons'].keys())} exons\n"
        f"  {len(primary_dat['locs'].keys())} seqids\n"
        f"  max_loc_length: {primary_dat['max_loc_length']}\n"
        f"  number filtered: {len(primary_dat['filtered'].keys())}\n\n"
    )


################################################
## if no secondary_gtf:

if not args.secondary_gtf:
    sys.stderr.write("No secondary_gtf.\n")

    sys.stderr.write(f"Generating final gtf output.\n")
    primary_list = output.primary_list(primary_dat)
    sys.stderr.write(f"  done at {round(time.time() - time_start, 2)} seconds\n\n")

    output_gtf = args.output_prefix + '.gtf'
    sys.stderr.write(f"Writing output_gtf {output_gtf}:\n")

    try:
        with open(output_gtf, 'w') as fh:
            for toks in primary_list:
                toks = [str(tok) for tok in toks]
                line = '\t'.join(toks)
                fh.write(f"{line}\n")
    except Exception as e:
        sys.stderr.write(f"ERROR: while writing output_gtf {output_gtf}:\n" f"{e}\n")
        sys.exit(35)

    sys.stderr.write(f"  done at {round(time.time() - time_start, 2)} seconds\n\n")
    sys.exit(0)


################################################
## parse secondary_gtf:

sys.stderr.write("Parsing secondary_gtf:\n")

try:
    with open(args.secondary_gtf, 'r') as fh:
        secondary_dat = gtf.parse_gtf(fh, args.secondary_prefix)

    sys.stderr.write(
        f"  done at {round(time.time() - time_start, 2)} seconds\n"
        f"  {len(secondary_dat['transcripts'].keys())} transcript_ids\n"
        f"  {len(secondary_dat['exons'].keys())} exons\n"
        f"  {len(secondary_dat['locs'].keys())} seqids\n"
        f"  max_loc_length: {secondary_dat['max_loc_length']}\n\n"
    )

except Exception as e:
    sys.stderr.write(
        f"ERROR: while processing secondary_gtf: '{args.secondary_gtf}':\n" 
        f"{e}\n"
    )
    sys.exit(37)


################################################
## resolve secondary_gtf against itself:

sys.stderr.write("Resolving secondary_gtf:\n")

secondary_dat2 = resolve.resolve_self(secondary_dat, 2, args)

sys.stderr.write(
    f"  done at {round(time.time() - time_start, 2)} seconds\n"
    f"  {len(secondary_dat2['transcripts'].keys())} transcript_ids\n"
    f"  {len(secondary_dat2['exons'].keys())} exons\n"
    f"  {len(secondary_dat2['locs'].keys())} seqids\n"
    f"  max_loc_length: {secondary_dat2['max_loc_length']}\n"
    f"  number filtered: {len(secondary_dat2['filtered'].keys())}\n\n"
)


################################################
## resolve resolved secondary against resolved primary:

sys.stderr.write("Resolving secondary_gtf against primary_gtf:\n")

try:
    output_dat = resolve.resolve_transcripts(primary_dat, secondary_dat2, args)
except Exception as e:
    sys.stderr.write(f"ERROR: while resolving:\n" f"{e}\n")
    sys.exit(43)

sys.stderr.write(
    f"  done at {round(time.time() - time_start, 2)} seconds\n"
    f"  {len(output_dat['transcripts'].keys())} transcript_ids\n"
    f"  {len(output_dat['exons'].keys())} exons\n"
    f"  {len(output_dat['locs'].keys())} seqids\n"
    f"  max_loc_length: {output_dat['max_loc_length']}\n"
    f"  number filtered: {len(output_dat['filtered'].keys())}\n\n"
)


#################################################################
## output:

sys.stderr.write(f"Generating final map and gtf outputs.\n")

mappings = output.filter_map(
    primary_dat['filtered'],
    secondary_dat2['filtered'],
    output_dat['filtered']
)

rmappings = output.reverse_map(
    mappings, 
    primary_dat['transcripts'], 
    secondary_dat2['transcripts']
)

## need to use original unresolved secondary_dat here:
output_dat = output.output_list(
    output_dat,
    rmappings,
    secondary_dat['transcripts'],
    args.ids_from
)
sys.stderr.write(f"  done at {round(time.time() - time_start, 2)} seconds\n\n")

output_tsv = args.output_prefix + ".filter.tsv"
sys.stderr.write(f"Writing output_tsv {output_tsv}:\n")

try:
    with open(output_tsv, 'w') as fh:
        labels = [ "filtered_db", "filtered_id", "kept_db", "kept_id" ]
        line = '\t'.join(labels)
        fh.write(f"{line}\n")
        for filtered_id in mappings:
            filtered_db, kept_db, kept_id = mappings[filtered_id]
            line = f"{filtered_db}\t{filtered_id}\t{kept_db}\t{kept_id}\n"
            fh.write(line)
except Exception as e:
    sys.stderr.write(f"ERROR: writing output_tsv {output_tsv}:\n" f"{e}\n")
    sys.exit(51)

output_gtf = args.output_prefix + ".gtf"
sys.stderr.write(f"Writing output_gtf {output_gtf}:\n")

try:
    with open(output_gtf, 'w') as fh:
        for toks in output_dat:
            toks = [str(tok) for tok in toks]
            line = '\t'.join(toks)
            fh.write(f"{line}\n")
except Exception as e:
    sys.stderr.write(f"ERROR: while writing output_gtf {output_gtf}:\n" f"{e}\n")
    sys.exit(53)

sys.stderr.write(f"  done at {round(time.time() - time_start, 2)} seconds\n\n")
sys.exit(0)

