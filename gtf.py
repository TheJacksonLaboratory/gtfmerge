#!/usr/bin/env python
# cython: embedsignature=True, binding=True

"""
parses gtf file; assumes input gtf has exon features that have 
  attributes transcript_id, and gene_id; 

main entrypoint is parse_gtf();

  takes filehandle and prefix to be prepended to gene_id and
    transcript_id.

  returns:
    locs: { seq: [ (start, end, strand, transcript), ...], ..., }
      locs[seq] sorted ascending by (start, end)
    transcripts: { transcript_id: [gene_id, seq, strand, start, end, [exon_keys]] }
      where exon_keys: 'seq:strand:start:end', and [exon_keys] sorted 
      ascending by (start, end)
    exons: { 'seq:strand:start:end': (seq, strand, start, end) }
    max_loc_length: length of longest encountered loc
"""

import re

######################################################################
## inferring features:

def fix_transcripts(exons, transcripts):

    """
    Inputs:

      exons: { 'seq:strand:start:end': (seq, strand, start, end) }

      transcripts: { 'transcript_id': [gene_id, seq, strand, start, end, [exon_keys]] }

    Output: None

    Side effects: for each transcripts[transcript_id]: sorts exons and fills start, end
    """

    def f_sort(exon_key):
        exon = exons[exon_key]
        start = exon[2]
        end = exon[3]
        return (start, end)

    for transcript_id in transcripts:
        exon_keys = transcripts[transcript_id][5]
        exon_keys.sort(key=f_sort)
        start_exon = exons[exon_keys[0]]
        end_exon = exons[exon_keys[-1]]
        transcripts[transcript_id][3] = start_exon[2]
        transcripts[transcript_id][4] = end_exon[3]


def infer_locs(transcripts):

    """
    Inputs:
      transcripts: { 'transcript_id': [gene_id, seq, strand, start, end, [exon_keys]] }

    Outputs:
      locs: { seq: [ (start, end, strand, transcript_id), ...], ..., }

      Note: locs[seq] are sorted by (start, end)
    """

    locs = {}

    for transcript_id in transcripts:
        transcript = transcripts[transcript_id]
        seq = transcript[1]
        if seq not in locs:
            locs[seq] = []

        ## (start, end, strand, transcript_id)
        locs[seq].append((transcript[3], transcript[4], transcript[2], transcript_id))

    for seq in locs:
        locs[seq].sort(key=lambda loc: (loc[0], loc[1]))

    return locs


def find_max_loc_length(locs):

    """
    Inputs:
      locs: { seq: [ (start, end, strand, transcript), ...], ..., }

    Output: length of longest loc; (1 + end - start)
    """

    max_loc_length = 0

    for seq in locs:
        for loc in locs[seq]:
            length = 1 + loc[1] - loc[0]
            if length > max_loc_length:
                max_loc_length = length

    return max_loc_length


######################################################################
## file parsing:

def register_exon(toks, exons, transcripts, attributes):

    """
    Description: Updates exons and transcripts for one exon record based on toks;
      exon defined by seq:strand:start:end; non-redundant collection saved.

    Inputs:
      toks: list of tokens, one per column from one row of input gtf file;
              0      1       2     3   4     5      6     7           8          9
        seqname source feature start end score strand frame [attributes] [comments]

        expects attributes to contain transcript_id and gene_id

      exons: { 'seq:strand:start:end': (seq, strand, start, end) }
        updated by this function;

      transcripts: { transcript_id: [gene_id, seq, strand, start, end, [exon_keys]] }
        updated by this function;

      attributes: dict w/ keys 'gene_id' and 'transcript_id'; values set to None if missing

    Output: None

    Side effect: updates exons and transcripts
    """

    if attributes['transcript_id'] is None:
        raise Exception(f"register_exon: could not find transcript_id in {toks}.")

    ## seq:strand:start:end:
    exon_key = f"{toks[0]}:{toks[6]}:{toks[3]}:{toks[4]}"

    exons[exon_key] = (
        toks[0],    ## 0:seq
        toks[6],    ## 1:strand
        toks[3],    ## 2:start
        toks[4],    ## 3:end
    )

    if attributes['transcript_id'] not in transcripts:
        transcripts[attributes['transcript_id']] = [
            attributes['gene_id'],                       ## 0: gene_id
            toks[0],                                     ## 1: seq
            toks[6],                                     ## 2: strand
            None,                                        ## 3: start
            None,                                        ## 4: end
            []                                           ## 5: exons
        ]

    transcripts[attributes['transcript_id']][5].append(exon_key)


def parse_attributes(tok9, prefix=None):

    """
    Description: parse gtf2 attributes (9th column) for gene_id and transcript_id.

    Inputs:

      tok9: 9th column (attributes) from one row from input gtf file.

      prefix: for transcript_id; can be None.

    Output: dict w/ keys 'gene_id' and 'transcript_id'; values set to None if missing,
    """

    if prefix is None:
        prefix = ""

    toks = re.split('; |;| ', tok9)
    toks = [ tok.strip('"') for tok in toks ]

    attributes = {}
    idx_last = len(toks) - 1
    idx = 0

    while(idx < idx_last):
        if toks[idx] == 'transcript_id':
            if 'transcript_id' in attributes:
                raise Exception(f"parse_attributes: multiple transcript_ids in {tok9}")
            elif idx == idx_last:
                raise Exception(f"parse_attributes: no value for transcript_id in {tok9}")
            else:
                attributes['transcript_id'] = f"{prefix}{toks[idx + 1]}"
                idx += 2
        elif toks[idx] == 'gene_id':
            if 'gene_id' in attributes:
                raise Exception(f"parse_attributes: multiple gene_ids in {tok9}")
            elif idx == idx_last:
                raise Exception(f"parse_attributes: no value for gene_id in {tok9}")
            else:
                attributes['gene_id'] = f"{prefix}{toks[idx + 1]}"
                idx += 2
        else:
            idx += 1

    return attributes


def parse_exon_line(line):

    """
    Input:
      line: one line from gtf file

    Output:
      toks is list of tokens, one per column from one row of input gtf file;
              0      1       2     3   4     5      6     7           8          9
        seqname source feature start end score strand frame [attributes] [comments]
    """

    toks = line.split('\t')
    toks = [ tok_i.strip() for tok_i in toks ]

    if len(toks) == 0:
        return None
    if toks[0].startswith('#'):
        return None
    if len(toks) < 9:
        raise Exception(f"parse_exon_line: fewer than 9 tokens on line.")

    ## start should be int:
    try:
        toks[3] = int(toks[3])
    except Exception as e:
        raise Exception(f"parse_exon_line: int() on '{toks[3]}': {e}")

    ## end should be int:
    try:
        toks[4] = int(toks[4])
    except Exception as e:
        raise Exception(f"parse_exon_line: int() on '{toks[4]}': {e}")

    return toks


def parse_exons(fh, prefix):

    """
    Inputs:

      fh: open readable filehandle to gtf;

      prefix: string to prepend to output transcript_ids;
        transcript_ids not changed if prefix is None.

    Outputs: exons dict:
      { 'seq:strand:start:end': (seq, strand, start, end) }
    """

    exons = {}
    transcripts = {}

    for line in fh:

        """
        toks is list of tokens, one per column from one row of input gtf file;
              0      1       2     3   4     5      6     7           8          9
        seqname source feature start end score strand frame [attributes] [comments]
        """

        try:
            toks = parse_exon_line(line)
        except Exception as e:
            raise Exception(f"parse_exons: {line}\n" f"parse_exon_line: {e}")

        if toks is None:
            continue

        if toks[2] != 'exon':
            continue

        try:
            attributes = parse_attributes(toks[8], prefix)
            register_exon(toks, exons, transcripts, attributes)
        except Exception as e:
            raise Exception(f"parse_exons: {line}\n" f"exon line error: {e}")

    return exons, transcripts


######################################################################
## main entry point:

def parse_gtf(fh, prefix):

    """
    Description: parses transcripts and exons from primary gtf;

    Inputs:

      fh: open readable filehandle to gtf;

      prefix: string to prepend to output transcript_ids;
        transcript_ids not changed if prefix is None.

    Outputs: dict w/ keys:
      locs: { seq: [ (start, end, strand, transcript), ...], ..., }
        locs[seq] sorted ascending by (start, end)
      transcripts: { transcript_id: [gene_id, seq, strand, start, end, [exon_keys]] }
        transcripts[transcript_id][5] (exon_keys) sorted ascending by (start, end)
      exons: { 'seq:strand:start:end': (seq, strand, start, end) }
      max_loc_length: length of longest encountered loc 
    """

    exons, transcripts = parse_exons(fh, prefix)
    fix_transcripts(exons, transcripts)
    locs = infer_locs(transcripts)
    max_loc_length = find_max_loc_length(locs)

    ## sort locs by start then end:
    for seq in locs:
        locs[seq].sort(key=lambda item: (item[0], item[1]))

    dat = {
        'locs': locs,
        'transcripts': transcripts,
        'exons': exons,
        'max_loc_length': max_loc_length
    }

    return dat

