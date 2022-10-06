#!/usr/bin/env python
# cython: embedsignature=True, binding=True

'''
used by mergegtfs.py

dat is dict w/ keys including:
  locs: { seq: [ [start, end, strand, transcript], ...], ..., }
  transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }
  exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }
  max_loc_length: length of longest transcript

params is a namespace w/ attributes including:
  tol_sj: tolerance for splice junction matching
  tol_tss: tolerance for start-site matching
  tol_tts: tolerance for termination-site matching
  nthreads: number of threads to use
  time_start: time (from time.time()) when job started
'''

def make_exon_recs(transcript, attributes, seq, dat):

    recs = []

    exon_keys = transcript[5]

    for exon_key in exon_keys:

        ## [ seq, strand, start, end]:
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

        recs.append(rec)

    return recs


def make_transcript_recs(transcript_id, seq, dat, ids):

    """
    ids['gene_ids']: { transcript_id1: gene_id1, transcript_id2: gene_id2, ... }
    ids['transcript_ids']: { transcript_id1: new_transcript_id1, transcript_id2: ... }
    """

    recs = []

    ##  [ gene_id, seq, strand, start, end, [exons] ]:
    transcript = dat['transcripts'][transcript_id]
    new_gene_id = ids['gene_ids'][transcript_id]
    new_transcript_id = ids['transcript_ids'][transcript_id]

    attributes = (
        f'gene_id "{new_gene_id}"; '
        f'transcript_id "{new_transcript_id}";'
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

    recs.append(rec)

    exon_recs = make_exon_recs(transcript, attributes, seq, dat)
    recs += exon_recs

    return recs


def dat2table(dat, ids):

    """
    ids['gene_ids']: { transcript_id1: gene_id1, transcript_id2: gene_id2, ... }
    ids['transcript_ids']: { transcript_id1: new_transcript_id1, transcript_id2: ... }
    """

    table = []       ## [ [seq,src,feature,start,end,score,strand,frame,attributes], ...]

    for seq in sorted(dat['locs']):

        ## [ [start, end, strand, transcript], ...]:
        locs = dat['locs'][seq]

        for loc in locs:
            transcript_id = loc[3]
            recs = make_transcript_recs(transcript_id, seq, dat, ids)
            table += recs

    return table


def xrefs(filtered, ids, gtfs):

    """
    Inputs:

      filtered: { filtered_id: kept_id, ... }

      ids: dict w/ two keys:
        ids['gene_ids']: { transcript_id1: gene_id1, transcript_id2: gene_id2, ... }
        ids['transcript_ids']: { transcript_id1: new_transcript_id1, transcript_id2: ... }

      gtfs: [input_gtf1, input_gtf2, ...]

    Output: [[source_gtf, filtered_id, kept_id, new_transcript_id, new_gene_id], ...]
    """

    table = []

    for old_id, kept_id in filtered.items():
        while kept_id in filtered:
            kept_id = filtered[kept_id]
        gtf_index, original_id = old_id.split('.', maxsplit=1)
        gtf_index = int(gtf_index)
        gtf_file = gtfs[gtf_index]
        gene_id = ids['gene_ids'][kept_id]
        transcript_id = ids['transcript_ids'][kept_id]
        table.append([gtf_file, original_id, transcript_id, gene_id])

    for old_id, gene_id in ids['gene_ids'].items():
        gtf_index, original_id = old_id.split('.', maxsplit=1)
        gtf_index = int(gtf_index)
        gtf_file = gtfs[gtf_index]
        transcript_id = ids['transcript_ids'][old_id]
        table.append([gtf_file, original_id, transcript_id, gene_id])

    table.sort(key=lambda it: (it[0], it[1]))
    return table


