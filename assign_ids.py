#!/usr/bin/env python
# cython: embedsignature=True, binding=True

'''
dat is dict w/ keys including:
  locs: { seq: [ (start, end, strand, transcript_id), ...], ..., }
  transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }
  exons: { 'seq:strand:start:end': ( seq, strand, start, end ) }
  max_loc_length: length of longest transcript

params is environment used to store gene_ids attribute.
'''

import output1

def exons_overlap(exon_keys1, exon_keys2, exons, params):

    n_exons1 = len(exon_keys1)
    n_exons2 = len(exon_keys2)
    last_idx = n_exons1 - 1
    n_match = 0

    for idx1 in range(n_exons1):

        ## [ seq, strand, start, end ]:
        exon1 = exons[exon_keys1[idx1]]
        len1 = 1 + exon1[3] - exon1[2]

        for idx2 in range(n_exons2):
            exon2 = exons[exon_keys2[idx2]]
            len2 = 1 + exon2[3] - exon2[2]
            length = min(len1, len2)
            if exon1[2] <= exon2[2]:
                overlap = exon1[3] - exon2[2]
            else:
                overlap = exon2[3] - exon1[2]
            if (overlap / length) >= params.p_exon_overlap:
                n_match += 1
            if (exon2[2] - exon1[3]) > 0:
                break

    if (n_match / min(n_exons1, n_exons2)) >= params.p_exons_overlap:
        return True
    else:
        return False


def match_genes(transcript1, transcript2, exons, params):

    if transcript1[2] != transcript2[2]:
        return False                       ## different strands

    exon_keys1 = transcript1[5]
    exon_keys2 = transcript2[5]

    return exons_overlap(exon_keys1, exon_keys2, exons, params)


def match_genes_after(idx, loc_list, dat, params, gene_ids):

    """
    adds to gene_ids
    """

    loc = loc_list[idx]          ## start, end, strand, transcript_id
    transcript_id = loc[3]
    transcript = dat['transcripts'][transcript_id]

    if transcript_id not in gene_ids:
        gene_ids[transcript_id] = f"{params.gene_prefix}{params.gene_idx}"
        params.gene_idx += 1

    for i in range(idx + 1, len(loc_list)):

        loc_i = loc_list[i]

        ## if start_i > end: no more overlaps
        if loc_i[0] > loc[1]:
            break

        ## if (start_i - end) > max_loc_length:
        ## if (loc_i[0] - loc[1]) > dat['max_loc_length']:
        ##     break                ## no more overlaps downstream

        transcript_id_i = loc_i[3]
        transcript_i = dat['transcripts'][transcript_id_i]

        if match_genes(transcript, transcript_i, dat['exons'], params):
            gene_ids[transcript_id_i] = gene_ids[transcript_id]


def assign_transcript_ids(loc_list, gene_ids, transcript_ids):

    """
    adds to transcript_ids
    loc_list: [ [start, end, strand, transcript_id], ...]
    gene_ids: { transcript_id1: gene_id1, transcript_id2: gene_id2, ... }
    transcript_ids: { transcript_id1: new_transcript_id1, transcript_id2: ... }
    assumes gene_ids unique to loc_list (chromosome); so may break w/ chrom fusions
    """

    indices = {}

    for loc  in loc_list:
        transcript_id = loc[3]
        gene_id = gene_ids[transcript_id]
        if gene_id not in indices:
            indices[gene_id] = 1
        else:
            indices[gene_id] += 1

        transcript_ids[transcript_id] = f"{gene_ids[transcript_id]}.{indices[gene_id]}"


def assign_ids(dat, params):

    """
    Inputs:
      dat: is dict with keys:
        locs: { seq: [ (start, end, strand, transcript_id), ...], ..., }
        transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }
        exons: { 'seq:strand:start:end': ( seq, strand, start, end ) }
        max_loc_length: length of longest transcript

      params: environment used to store attribute params.gene_idx, used 
        to index successive genes.

    Output: dict w/ keys:
      'gene_ids': { transcript_id1: gene_id1, transcript_id2: gene_id2, ... }
      'transcript_ids': { transcript_id1: new_transcript_id1, transcript_id2: ...}
    """

    gene_ids = {}
    transcript_ids = {}
    params.gene_idx = 1

    for seq in sorted(dat['locs']):

        print(f"assign_ids: start seq {seq} at {output1.elapsed(params)} seconds")

        loc_list = dat['locs'][seq]

        for idx in range(len(loc_list)):
            match_genes_after(idx, loc_list, dat, params, gene_ids)

        assign_transcript_ids(loc_list, gene_ids, transcript_ids)

    ids = {
      'gene_ids': gene_ids,
      'transcript_ids': transcript_ids
    }

    return ids

