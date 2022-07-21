'''
dat is dict w/ keys including:
  locs: { seq: [ [start, end, strand, transcript_id], ...], ..., }
  transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }
  exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }
  max_loc_length: length of longest transcript

params is dict w/ keys including:
  tol_sj: tolerance for splice junction matching
  tol_tss: tolerance for start-site matching
  tol_tts: tolerance for termination-site matching
  nthreads: number of threads to use
  time_start: time (from time.time()) when job started
'''

def exons_overlap(exon_keys1, exon_keys2, exons, cutoff=0.5):

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
            if (overlap / length) >= cutoff:
                n_match += 1
            if (exon2[2] - exon1[3]) > 0:
                break

    if (n_match / min(n_exons1, n_exons2)) >= cutoff:
        return True
    else:
        return False


def match_genes(transcript1, transcript2, exons, params):

    if transcript1[2] != transcript2[2]:
        return False                       ## different strands

    strand = transcript1[2]

    if strand == '+':
        tol_first = params.tol_tss
        tol_last = params.tol_tts
    elif strand == '-':
        tol_first = params.tol_tts
        tol_last = params.tol_tss
    else:
        raise Exception(
            f"get_tolerances: unexpected strand '{strand}'"
        )

    tol_sj = params.tol_sj

    exon_keys1 = transcript1[5]
    exon_keys2 = transcript2[5]

    return exons_overlap(exon_keys1, exon_keys2, exons, cutoff=0.5)


def match_genes_after(idx, loc_list, dat, params, gene_ids):

    loc = loc_list[idx]          ## start, end, strand, transcript_id
    transcript_id = loc[3]
    transcript = dat['transcripts'][transcript_id]

    if transcript_id not in gene_ids:
        gene_ids[transcript_id] = f"PB.{params.gene_idx}"
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


def resolve_ids(dat, params):

    """
    ids['gene_ids']: { transcript_id1: gene_id1, transcript_id2: gene_id2, ... }
    ids['transcript_ids']: { transcript_id1: new_transcript_id1, transcript_id2: ...}
    """

    gene_ids = {}
    transcript_ids = {}
    params.gene_idx = 1

    for seq in dat['locs']:

        print(f"DEBUG: start seq: {seq} at {output.elapsed(params)} seconds")

        loc_list = dat['locs'][seq]

        for idx in range(len(loc_list)):
            match_genes_after(idx, loc_list, dat, params, gene_ids)
        print(f"DEBUG: match_genes_after done at {output.elapsed(params)} seconds")

        assign_transcript_ids(loc_list, gene_ids, transcript_ids)
        print(f"DEBUG: assign_transcript_ids done at {output.elapsed(params)} seconds")

    ids = {
      'gene_ids': gene_ids,
      'transcript_ids': transcript_ids
    }

    return ids

