#!/usr/bin/env python
# cython: embedsignature=True, binding=True

"""
Code for identifying overlaps between transcripts and exons.
"""

import bisect


######################################################################################
## matching exons and transcripts:

def exons_match(exon1, exon2, tol_start, tol_end):

    '''
    Description: compares exons of two transcripts, returns True if match else False.

    Inputs:
      exon1: ( seq, strand, start, end )
      exon2: ( seq, strand, start, end )
      tol_start: tolerance for matching exon start coordinates.
      tol_end: tolerance for matching exon end coordinates.

    Output: returns True if match within tolerances; otherwise returns False.

    NOTE: start and end are relative to sequence start on positive strand;
      so need to flip tolerances for negative strand matching of terminal exons.
    '''

    if exon1[0] != exon2[0]:
        return False                             ## seqs don't match

    if exon1[1] != exon2[1]:
        return False                             ## strand doesn't match

    dist = abs(exon1[2] - exon2[2])

    if dist > tol_start:
        return False                             ## starts don't match

    dist = abs(exon1[3] - exon2[3])
    if dist > tol_end:
        return False                             ## ends don't match

    return True                                  ## everything w/i tolerances


def match_transcripts(transcript1, transcript2, exons1, exons2, args):

    """
    Description: 

    Inputs:
      transcript1: [ gene_id, seq, strand, start, end, [exons] ]

      transcript2: [ gene_id, seq, strand, start, end, [exons] ]

      exons1: { 'seq:strand:start:end': [ seq, strand, start, end ], ... }
        corresponding to transcript1[5] exon_keys

      exons2: { 'seq:strand:start:end': [ seq, strand, start, end ], ... }
        corresponding to transcript2[5] exon_keys

      args is dict w/ keys including:
        tol_sj: tolerance for splice junction matching;
        tol_tss: tolerance for start-site matching;
        tol_tts: tolerance for termination-site matching;

    Returns True if transcripts match, else False.
    """

    if transcript1[2] != transcript2[2]:
        return False                       ## different strands

    if len(transcript1[5]) != len(transcript2[5]):
        return False                       ## different number of exons

    if transcript1[2] == '+':
        tol_first = args.tol_tss
        tol_last = args.tol_tts
    elif transcript1[2] == '-':
        tol_first = args.tol_tts
        tol_last = args.tol_tss
    else:
        raise Exception(
            f"match_transcripts: unexpected strand '{transcript1[2]}' "
            f"for transcript {transcript1}."
        )

    tol_sj = args.tol_sj

    ## iterate over exons:

    exon_keys1 = transcript1[5]
    exon_keys2 = transcript2[5]
    n_exons = len(exon_keys1)
    last_idx = n_exons - 1

    for idx in range(n_exons):
        if(n_exons == 1):
            tol_start = tol_first
            tol_end = tol_last
        elif(idx == 0):
            tol_start = tol_first
            tol_end = tol_sj
        elif(idx == last_idx):
            tol_start = tol_sj
            tol_end = tol_last
        else:
            tol_start = tol_sj
            tol_end = tol_sj

        exon1 = exons1[exon_keys1[idx]]
        exon2 = exons2[exon_keys2[idx]]

        if not exons_match(exon1, exon2, tol_start, tol_end):
            return False

    return True


######################################################################################
## transfer data:

def update_entries(seq, filtered, source, dat1, dat2):

    """
    Description: updates dat2 w/ info in dat1 corresponding to locs in 
      dat1['locs'][seq] that are not filtered;

    Inputs:
      seq: seq key for dat1['locs'];

      filtered: info on filtered locs in dat1['locs'][seq]:
        { filtered_id: kept_id, ... }; 

      source: 1 if dat1 represents primary_gtf, 2 if represents secondary_gtf.

      dat1 is dict w/ keys:
        locs: { seq: [ [start, end, strand, transcript], ...], ..., }
        transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }
        exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }
        max_loc_length: length of longest transcript

      dat2: dict w/ keys like dat1:
        locs, transcripts, exons, and max_loc_length; 
        also, filtered transcript_ids in filtered dict: { filtered_id: kept_id, ... }
        and sources of kept transcript_ids in sources: { kept_id: source, ... }
          where source is 1 for primary_gtf, 2 for secondary_gtf;

    Output: None;

    Side effects: 
      updates dat2 w/ info from dat1 for locs in dat1['locs'][seq] that are not filtered;
        ensures locs stays sorted;
      updates dat2['filtered'] w/ filtered transcript_ids from loc_list;
      updates dat2['sources'] for locs in dat1['locs'][seq] that are not filtered;
        { kept_id: source, ... }; error if kept_id duplicated.
    """

    if seq not in dat2['locs']:
        dat2['locs'][seq] = []

    loc_list = dat1['locs'].get(seq)    ## [ [start, end, strand, transcript], ...]
    if loc_list is None:
        return

    for loc in loc_list:    ## iterate over loc_list from dat1

        transcript_id = loc[3]

        if transcript_id in filtered:
            dat2['filtered'][transcript_id] = filtered[transcript_id]
            continue 

        if transcript_id in dat2['sources']:
            raise Exception(f"update_entries: already saw transcript_id '{transcript_id}'")

        dat2['sources'][transcript_id] = source
        dat2['locs'][seq].append(loc)

        transcript = dat1['transcripts'][transcript_id]
        dat2['transcripts'][transcript_id] = transcript

        exon_keys = transcript[5]
        for exon_key in exon_keys:
            dat2['exons'][exon_key] = dat1['exons'][exon_key]

        ## start - end:
        length = 1 + loc[1] - loc[0]

        if length > dat2['max_loc_length']:
            dat2['max_loc_length'] = length

    ## ensure locs stay sorted:
    dat2['locs'][seq].sort(key=lambda item: (item[0], item[1]))


######################################################################################
## filter against self:

def filter_after(idx, loc_list, dat, args):

    """
    Description:

    Usage: filtered = filter_after(idx, loc_list, args)

    Inputs:
      idx: integer index in [0, Inf) of current loc in loc_list

      loc_list: [ [start, end, strand, transcript], ...]

      dat is dict w/ keys:
        locs: { seq: [ [start, end, strand, transcript], ...], ..., }
        transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }
        exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }
        max_loc_length: length of longest transcript

      args is dict w/ keys including:
        tol_sj: tolerance for splice junction matching;
        tol_tss: tolerance for start-site matching;
        tol_tts: tolerance for termination-site matching;

    Outputs:
      dict w/ info on locs in loc_list to be filtered:
        { filtered_id: kept_id, ... }
      dict can be empty if none.
    """

    filtered = {}
    loc = loc_list[idx]          ## start, end, strand, transcript_id
    transcript_id = loc[3]
    transcript = dat['transcripts'][transcript_id]

    for i in range(idx + 1, len(loc_list)):

        loc_i = loc_list[i]

        ## if start_i > end: no more matches:
        if loc_i[0] > loc[1]:
            break

        ## if (start_i - end) > max_loc_length:
        ## if (loc_i[0] - loc[1]) > dat['max_loc_length']:
        ##     break                ## no more overlaps downstream

        transcript_id_i = loc_i[3]
        transcript_i = dat['transcripts'][transcript_id_i]

        if match_transcripts(transcript, transcript_i, dat['exons'], dat['exons'], args):
            if transcript_id_i not in filtered:
                filtered[transcript_id_i] = transcript_id

    return filtered


def resolve_self(dat, source, args):

    """
    Description: collapses redundant transcripts.

    Inputs:

      dat is dict w/ keys:
        locs: { seq: [ [start, end, strand, transcript], ...], ..., }
        transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }
        exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }
        max_loc_length: length of longest transcript

      source: 1 if dat from primary_gtf, 2 if from secondary_gtf;

      args is dict w/ keys including:
        tol_sj: tolerance for splice junction matching;
        tol_tss: tolerance for start-site matching;
        tol_tts: tolerance for termination-site matching;

    Returns dict w/ info on non-redundant transcripts in keys:
      locs, transcripts, exons, and max_loc_length; 
      also, filtered: { filtered_id: kept_id, ... }
      and sources: { kept_id: source, ... }
        where source is 1 for primary_gtf, 2 for secondary_gtf;
      Redundancy determined using tol_sj, tol_tss, and tol_tts;
    """

    output = {
        'locs': {},
        'transcripts': {},
        'exons': {},
        'max_loc_length': 0,
        'filtered': {},
        'sources': {}
    }

    for seq in dat['locs']:

        loc_list = dat['locs'][seq]
        filtered = {}             ## { filtered_id: kept_id, ... }

        for idx in range(len(loc_list)):
            ## { filtered_id: kept_id, ... }
            filtered_i = filter_after(idx, loc_list, dat, args) 
            for transcript_id in filtered_i:
                if transcript_id not in filtered:
                    filtered[transcript_id] = filtered_i[transcript_id]

        update_entries(seq, filtered, source, dat, output)

    return output


######################################################################################
## filter one set against another:

def filter_locs(seq, primary_dat, secondary_dat, args):

    """
    Inputs:

      seq: key into *_dat['locs']

      primary_dat is a dict w/ keys including:
        locs: { seq: [ [start, end, strand, transcript], ...], ..., }
        transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }
        exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }
        max_loc_length: length of longest transcript

      secondary_dat is a dict w/ keys like primary_dat:
        locs, transcripts, exons, max_loc_length.

      args is dict w/ keys including:
        tol_sj: tolerance for splice junction matching;
        tol_tss: tolerance for start-site matching;
        tol_tts: tolerance for termination-site matching;

    Outputs:
      dict mapping filtered transcripts to unfiltered transcript:
        { filtered_id: kept_id, ... }
    """

    filtered = {}

    ## [ [start, end, strand, transcript], ...]:
    primary_locs = primary_dat['locs'].get(seq)
    secondary_locs = secondary_dat['locs'].get(seq)

    if primary_locs is None:
        return filtered

    if secondary_locs is None:
        return filtered

    max_loc_length = max(
        primary_dat['max_loc_length'], 
        secondary_dat['max_loc_length']
    )

    ## for bisect search; accommodates older python w/o bisect_left 'key' parameter:
    secondary_starts = [loc[0] for loc in secondary_locs]

    ##  [start, end, strand, transcript]:

    for primary_loc in primary_locs:

        primary_id = primary_loc[3]
        primary_transcript = primary_dat['transcripts'][primary_id]
        left_pos = primary_loc[0] - max_loc_length

        if left_pos <= 0:
            left_pos = 0
            idx = 0
        else:
            ## assumes secondary_locs sorted by start position:
            idx = bisect.bisect_left(secondary_starts, left_pos)
            if idx == len(secondary_locs):
                continue         ## no potential matches; done w/ primary_loc

        for secondary_loc in secondary_locs[idx:]:

            ## secondary_start - primary_start; after potential matches:
            if (secondary_loc[0] - primary_loc[0]) > max_loc_length:
               break             ## after potential matches; done w/ primary_loc

            secondary_id = secondary_loc[3]
            secondary_transcript = secondary_dat['transcripts'][secondary_id]

            if match_transcripts(
                primary_transcript, 
                secondary_transcript, 
                primary_dat['exons'],
                secondary_dat['exons'],
                args
            ):
                filtered[secondary_id] = primary_id

    return filtered


def resolve_transcripts(primary_dat, secondary_dat, args):

    """
    Inputs:

      primary_dat is a dict w/ keys including:
        locs: { seq: [ [start, end, strand, transcript], ...], ..., }
        transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }
        exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }
        max_loc_length: length of longest transcript

      secondary_dat is a dict w/ keys like primary_dat: 
        locs, transcripts, exons, max_loc_length.

      args is dict w/ keys including:
        tol_sj: tolerance for splice junction matching;
        tol_tss: tolerance for start-site matching;
        tol_tts: tolerance for termination-site matching;

    Output: returns dict w/ info on non-redundant transcripts
      in secondary dat w/ keys like primary_dat: 
        locs, transcripts, exons, max_loc_length
      but also filtered w/ info on filtered transcripts:
        { filtered_id: kept_id, ... }
      and sources w/ info on source of kept transcripts:
        { kept_id: source, ...}; where source is 1 for primary_gtf, else 2.
    """

    output_dat = {
        'locs': {},
        'transcripts': {},
        'exons': {},
        'max_loc_length': 0,
        'filtered': {},
        'sources': {}
    }

    for seq in primary_dat['locs']:

        ## keep all entries from primary_dat (nothing filtered):
        update_entries(seq, [], 1, primary_dat, output_dat)

        ## find redundant entries in secondary_dat:
        ##   { filtered_id: kept_id, ... }
        filtered = filter_locs(seq, primary_dat, secondary_dat, args)

        ## keep only non-redundant entries from secondary_dat:
        update_entries(seq, filtered, 2, secondary_dat, output_dat)

    ## supplement w/ seqs unique to secondary_dat:
    for seq in secondary_dat['locs']:
        if seq in primary_dat['locs']:
            continue
        update_entries(seq, [], 2, secondary_dat, output_dat)

    ## sort locs by start then end:
    for seq in output_dat['locs']:
        output_dat['locs'][seq].sort(key=lambda item: (item[0], item[1]))

    return output_dat

