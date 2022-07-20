#!/usr/bin/env python3

## system:
import multiprocessing
import sys
import time

## local:
import initial2
import gtf
import output
import resolve

'''
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

#########################################################################
## output:


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


##########################################################################
## assiging gene_ids + transcript_ids:

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


##########################################################################
## resolving transcripts:

def read_gtf_list_file(gtf_list_file):
    '''
    inputs:
        gtf_list_file: a text file w/ one path to gtf file per line

    output:
        list of gtf file paths
    '''

    gtf_files = []

    try:
        with open(gtf_list_file, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not len(line):
                    continue
                if line.startswith('#'):
                    continue
                gtf_files.append(line)
    except Exception as e:
        raise Exception(f"read_gtf_list_file: for gtf_list_file {gtf_list_file}: {e}")

    return gtf_files


def load_gtf(gtf_path, prefix, params):

    try:
        with open(gtf_path, 'r') as fh:
            dat = gtf.parse_gtf(fh, prefix)
    except Exception as e:
        raise Exception(f"load_gtf: processing {gtf_path}: {e}")

    dat = resolve.resolve_self(dat, 1, params)
    return dat


def pop_chunk(setup, chunk_size):

    '''
    changes setup
    '''

    chunk = []
    length = len(setup)

    if length:
        idx_stop = chunk_size
        if idx_stop > length:
            idx_stop = length
        chunk += setup[0 : idx_stop]
        del setup[0 : idx_stop]

    return chunk


def collapse_list(dat_list, params, pool):

    n = len(dat_list)

    n_pairs = int(n / 2)
    last_paired = True
    if n % 2:
        last_paired = False
        n_pairs += 1
    idx_last = n_pairs - 1

    setup = []

    for idx in range(n_pairs):
        i = idx * 2
        if idx == idx_last:
            if last_paired:
                setup.append([dat_list[i], dat_list[i+1], params])
            else:
                pass
        else:
           setup.append([dat_list[i], dat_list[i+1], params])

    print(f"len(setup): {len(setup)}")

    ## resolve_transcripts(primary_dat, secondary_dat, params):
    if setup:
        dat_collapse = pool.starmap(resolve.resolve_transcripts, setup)
    else:
        dat_collapse = []

    if not last_paired:
        dat_collapse.append(dat_list[idx_last * 2])

    return dat_collapse


def process_gtfs(params, pool):

    gtfs = read_gtf_list_file(params.gtf_list_file)
    setup = [(gtf, f"{i}.", params) for (gtf, i) in zip(gtfs, range(len(gtfs)))]

    dat = None

    chunk = pop_chunk(setup, params.nthreads)
    while chunk:
        print(f"len(chunk): {len(chunk)}")
        dat_list = pool.starmap(load_gtf, chunk)
        print(f"process_gtfs: load_gtf done at {output.elapsed(params)} seconds")
        print(f"len(dat_list): {len(dat_list)}")
        while len(dat_list) > 1:
            dat_list = collapse_list(dat_list, params, pool)
            print(
                f"process_gtfs: collapse_list done at "
                f"{output.elapsed(params)} seconds"
            )
            print(f"len(dat_list): {len(dat_list)}")

        if dat:
            dat = resolve.resolve_transcripts(dat, dat_list[0], params)    
            print(
                f"process_gtfs: resolve_transcripts done at "
                f"{output.elapsed(params)} seconds"
            )
        else:
            dat = dat_list[0]

        chunk = pop_chunk(setup, params.nthreads)

    return dat


##########################################################################
## main:

if __name__ == '__main__':

    params = initial2.initialize()

    try:
        with multiprocessing.Pool(params.nthreads) as pool:
            dat = process_gtfs(params, pool)
    except Exception as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(31)

    print(f"main: process_gtfs done at {output.elapsed(params)} seconds")

    ids = resolve_ids(dat, params)
    print(f"main: resolve_ids done at {output.elapsed(params)} seconds")

    gtf_table = dat2table(dat, ids)
    print(f"main: dat2table done at {output.elapsed(params)} seconds")

    output_gtf = f"{params.output_prefix}.gtf"
    try:
        with open(output_gtf, 'w') as fh:
            for toks in gtf_table:
                toks = [str(tok) for tok in toks]
                line = '\t'.join(toks)
                fh.write(f"{line}\n")
    except Exception as e:
        sys.stderr.write(f"ERROR: writing {output_gtf}: {e}\n")
        sys.exit(33)

    print(
        f"main: finished writing {output_gtf} at "
        f"{output.elapsed(params)} seconds"
    )

    sys.exit(0)


