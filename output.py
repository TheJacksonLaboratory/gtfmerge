#!/usr/bin/env python
# cython: embedsignature=True, binding=True

"""
Code for generating and printing outputs;

used by gtfmerge.py
"""

import copy
import re
import time

def elapsed(params, places=3):
    return round(time.time() - params.time_start, places)


def time_stamp():

    """
    returns human readable sortable timestamp str
    """

    time_stamp = time.strftime(
        '%Y%m%d%H%M%S',
        time.localtime(time.time())
    )

    return time_stamp


def filter_map(primary_filter, secondary_filter, output_filter):

    """
    Inputs:
      primary_filter: {filtered_primary_id: kept_primary_id, ...}
      secondary_filter: {filtered_secondary_id: kept_secondary_id, ... }
      output_filter: {filtered_secondary_id: kept_primary_id, ... }

    Output: { filter_id: [filter_db, kept_db, kept_id], ... }
    """

    mappings = {}    ## { filter_id: [filter_db, kept_db, kept_id], ... }

    for filtered_id in primary_filter:

        kept_id = primary_filter[filtered_id]
        while kept_id in primary_filter:
            kept_id = primary_filter[kept_id]

        ## sourcedb, xrefdb, xref_id
        mappings[filtered_id] = [1, 1, kept_id]

    for filtered_id in secondary_filter:

        kept_id = secondary_filter[filtered_id]
        while kept_id in primary_filter:
            kept_id = secondary_filter[kept_id]

        mappings[filtered_id] = [2, 2, kept_id]

    for filtered_id in output_filter:

        kept_id = output_filter[filtered_id]
        mappings[filtered_id] = [2, 1, kept_id]

    return mappings


def reverse_map(mappings, primary_transcripts, secondary_transcripts):

    """
    Description: generates mappings from kept transcript_ids to 
        filtered transcript_ids

    Inputs:

      mappings: { filter_id: [filter_db, kept_db, kept_id], ... }
        kept_db in {1, 2}, where 2 is secondary_gtf

      *_transcripts: info on transcripts:
        { 'transcript_id': [ gene_id, seq, strand, start, end, [exons]] }
    
    Output: mappings from kept transcripts to list of filtered
      { kept_id: [ [kept_db, filtered_db, filtered_id], ... ], ... }
    """

    rev_map = {}

    for filtered_id in mappings:

        filtered_db, kept_db, kept_id = mappings[filtered_id]

        if kept_id not in rev_map:
            rev_map[kept_id] = []

        rev_map[kept_id].append([kept_db, filtered_db, filtered_id])

    return rev_map


def best_secondary(kept_id, rev_map, secondary_transcripts):

    """
    Inputs:

      kept_id: transcript_id for kept transcript;

      rev_map: mappings from kept transcripts to list of filtered
        ones from secondary_gtf, ordered by length (1 + end - start):
        { kept_id: [ [kept_db, filtered_db, filtered_id], ... ], ... }

      secondary_transcripts: input data from secondary_gtf; dict:
        { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }

    Output: transcript_id of longest filtered secondary transcript
      corresponding to kept_id; or empty str '' if no corresponding 
      secondary transcript;
    """

    id_list = rev_map.get(kept_id)

    if not id_list:
        return ""

    ## keep only filtered_id where filtered_db == 2:
    id_list = [item[2] for item in id_list if item[1] == 2]

    if not id_list:
        return ""
    
    def transcript_length(transcript_id):
        transcript = secondary_transcripts[transcript_id]
        length = 1 + transcript[4] - transcript[3]
        return length

    id_list.sort(key=transcript_length, reverse=True)
    return id_list[0]


def make_attributes(output_id, output_dat, args, ids=None, secondary_transcripts=None, rev_map=None, ids_from=1):

    """
    Inputs:
      output_id: a kept transcript_id

      output_dat: output data; dict w/ keys including:
        locs: { seq: [ [start, end, strand, transcript], ...], ..., }
        transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ], ... }
        exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }
        max_loc_length: length of longest transcript
        sources: { kept_id: source, ... }; where source is 1 for primary_gtf, else 2;

      args: environment with attributes including primary_prefix and secondary_prefix

      ids: dict w/ keys:
        'gene_ids': { transcript_id1: gene_id1, transcript_id2: gene_id2, ... }
        'transcript_ids': { transcript_id1: new_transcript_id1, transcript_id2: ...}
        if None, no gene remapping.

      secondary_transcripts: input data from secondary_gtf; dict:
        { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ], ... }
        if None, no secondary transcript_ids etc.

      rev_map: mappings from kept transcripts to list of filtered
        ones from secondary_gtf, ordered by length (1 + end - start):
          { kept_id: [ [kept_db, filtered_db, filtered_id], ... ], ... }
        if None, no secondary transcript_ids etc.

      ids_from: where to get gene_id and transcript_id from; 1: primary_gtf; 2: secondary;

    Output: attributes string with gene_id, transcript_id, primary_id and secondary_id:
      'gene_id "gene328"; transcript_id "rna5783.2"; primary_id "rna5783.2"; secondary_id "PB.1.2";'
    """

    ## [ gene_id, seq, strand, start, end, [exon_keys] ]
    output_transcript = output_dat['transcripts'][output_id]
    source = output_dat['sources'][output_id]

    if source == 1:
        primary_id = output_id
        if rev_map and secondary_transcripts:
            secondary_id = best_secondary(output_id, rev_map, secondary_transcripts)
        else:
            secondary_id = ''
    else:
        primary_id = ""
        secondary_id = output_id

    if secondary_id and ids_from == 2:
        secondary_transcript = secondary_transcripts[secondary_id]
        gene_id = secondary_transcript[0]
        transcript_id = secondary_id
    else:
        gene_id = output_transcript[0]
        transcript_id = output_id

    if args.primary_prefix:
        primary_id = primary_id.lstrip(args.primary_prefix)

    if args.secondary_prefix:
        secondary_id = secondary_id.lstrip(args.secondary_prefix)

    if ids:
        gene_id = ids['gene_ids'][output_id]

    output = (
        f'gene_id "{gene_id}"; '
        f'transcript_id "{transcript_id}"; '
        f'primary_id "{primary_id}"; '
        f'secondary_id "{secondary_id}";'
    )

    return output


def primary_list(dat, args):

    """
    Inputs:

      dat is set of transcripts to save; a dict w/ keys including:
        locs: { seq: [ [start, end, strand, transcript], ...], ..., }
        transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }
        exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }
        max_loc_length: length of longest transcript
        sources: { kept_id: source, ... }; where source is 1 for primary_gtf, else 2;

      args: environment with attributes including primary_prefix and secondary_prefix.

    Output: list w/ format:
      [seq source feature start end score strand frame attributes]
      where attributes has gene_id, transcript_id, primary_id and secondary_id:
        gene_id "gene328"; transcript_id "rna5783.2"; primary_id "rna5783.2"; secondary_id ""
    """

    output = []

    for seq in dat['locs']:

        ## [ [start, end, strand, transcript], ...]:
        locs = dat['locs'][seq]

        for loc in locs:
            transcript_id = loc[3]

            ##  [ gene_id, seq, strand, start, end, [exons] ]:
            transcript = dat['transcripts'][transcript_id]

            ## secondary_transcripts and rev_map None; source = 1:
            attributes = make_attributes(transcript_id, dat, args, None, None, None, 1)

            rec = [
                seq,                ## sequence/chromosome
                'gtfmerge',         ## source
                'transcript',       ## feature
                transcript[3],      ## start
                transcript[4],      ## end
                '.',                ## score
                transcript[2],      ## strand
                '.',                ## frame
                attributes          ## attributes
            ]

            output.append(copy.deepcopy(rec))

            exon_keys = transcript[5]
            for exon_key in exon_keys:
                ## [ seq, strand, start, end ]
                exon = dat['exons'][exon_key]

                rec = [
                    seq,            ## sequence/chromosome
                    'gtfmerge',     ## source
                    'exon',         ## feature
                    exon[2],        ## start
                    exon[3],        ## end
                    '.',            ## score
                    exon[1],        ## strand
                    '.',            ## frame
                    attributes      ## attributes
                ]

                output.append(copy.deepcopy(rec))

    return output


def output_list(dat, rev_map, secondary_transcripts, ids, args):

    """
    Inputs:

      dat is set of transcripts to save; a dict w/ keys including:
        locs: { seq: [ [start, end, strand, transcript], ...], ..., }
        transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }
        exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }
        max_loc_length: length of longest transcript
        sources: { kept_id: source, ... }; where source is 1 for primary_gtf, else 2;

      rev_map: mappings from kept transcripts to list of filtered
        ones from secondary_gtf, ordered by length (1 + end - start):
        { kept_id: [ filter_id, ... ], ... }

      secondary_transcripts:  
        { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }

      ids: dict w/ keys:
        'gene_ids': { transcript_id1: gene_id1, transcript_id2: gene_id2, ... }
        'transcript_ids': { transcript_id1: new_transcript_id1, transcript_id2: ...}

      args: environment with attributes including ids_from, primary_prefix and secondary_prefix
        ids_from: in {1, 2}; if 1, get gene_id and transcript_id from primary_gtf,
          else from secondary_gtf.

    Output: list:
      [seq source feature start end score strand frame attributes]
      where attributes has gene_id, transcript_id, primary_id and secondary_id:
        gene_id "gene328"; transcript_id "rna5783.2"; primary_id "rna5783.2"; secondary_id "PB.1.2"  
    """

    output = []       ## [ [seq,src,feature,start,end,score,strand,frame,attributes], ...]

    for seq in sorted(dat['locs']):

        ## [ [start, end, strand, transcript], ...]:
        locs = dat['locs'][seq]

        for loc in locs:
            transcript_id = loc[3]

            ##  [ gene_id, seq, strand, start, end, [exons] ]:
            transcript = dat['transcripts'][transcript_id]

            attributes = make_attributes(
                transcript_id, 
                dat, 
                args,
                ids,
                secondary_transcripts, 
                rev_map, 
                args.ids_from,
            ) 

            rec = [
                seq,                ## sequence/chromosome
                'gtfmerge',         ## source
                'transcript',       ## feature
                transcript[3],      ## start
                transcript[4],      ## end
                '.',                ## score
                transcript[2],      ## strand
                '.',                ## frame
                attributes          ## attributes
            ]

            output.append(copy.deepcopy(rec))

            exon_keys = transcript[5]
            for exon_key in exon_keys:
                ## [ seq, strand, start, end ]
                exon = dat['exons'][exon_key]

                rec = [
                    seq,            ## sequence/chromosome
                    'gtfmerge',     ## source
                    'exon',         ## feature
                    exon[2],        ## start
                    exon[3],        ## end
                    '.',            ## score
                    exon[1],        ## strand
                    '.',            ## frame
                    attributes      ## attributes
                ]

                output.append(copy.deepcopy(rec))
                
    return output 


def write_output(output_list, out_file):

    """
    Inputs:

      output_list is a list of row data:
        [ [seq, source, feature, start, end, score, strand, frame, attributes], ... ]

      out_file: path to file where output should be written;
        overwrites if exists;

    Output: None

    Side effects: prints tab-delited output with column format:

              0      1       2     3   4     5      6     7            8 
        seqname source feature start end score strand frame [attributes]

        where attributes has gene_id, transcript_id, primary_id and secondary_id: 
          gene_id "gene328"; transcript_id "rna5783.2"; primary_id "rna5783.2"; secondary_id "PB.1.2"
    """

    try:
        with open(out_file, 'w') as fh:
            for toks in output_list:
                line = '\t'.join(toks)
                fh.write(f"{line}\n")
    except Exception as e:
        raise Exception(f"write_output: while writing out_file '{out_file}': {e}")


def write_map(mappings, map_file):

    """
    Inputs:
      mappings: { filter_id: [filter_db, kept_db, kept_id], ... }

      map_file: path to file where mappings should be written; 
        overwrites if exists;

    Output: None

    Side effects: prints final mappings to map_file; line format:
      f"{filter_db}\t{filter_id}\t{kept_db}\t{kept_id}\n"
    """

    try:
        with open(map_file, 'w') as fh:
            labels = [ "filtered_db", "filtered_id", "kept_db", "kept_id" ]
            line = '\t'.join(labels)
            fh.write(f"{line}\n")
            for filtered_id in mappings:
                filtered_db, kept_db, kept_id = mappings[filtered_id]
                line = f"{filtered_db}\t{filtered_id}\t{kept_db}\t{kept_id}\n"
                fh.write(line)
    except Exception as e:
        raise Exception("write_map: while writing to map_file '{map_file}': {e}")

