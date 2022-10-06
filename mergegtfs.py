#!/usr/bin/env python3

## system:
import ray
import sys
import time

## local:
import assign_ids
import initial2
import gtf
import output1
import output2
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

    gtf_files = set(gtf_files)
    gtf_files = list(gtf_files)
    gtf_files.sort()

    return gtf_files


@ray.remote
def load_gtf_ray(gtf_path, prefix, params):

    try:
        with open(gtf_path, 'r') as fh:
            dat = gtf.parse_gtf(fh, prefix)
    except Exception as e:
        raise Exception(f"load_gtf: processing {gtf_path}: {e}")

    dat = resolve.resolve_self(dat, 1, params)
    return dat


def pop_chunk(setup, chunk_size):

    '''
    inputs: 
      setup:  [(input_gtf1, 1), (input_gtf2, 2), ..., (input_gtfN, N)]
      chunk_size: number of records per chunk

    ouput: chunk of requested size or consisting of remaining records,
      whichever is less; [(gtf1, 1), (gtf2, 2), ... ]

    side-effects: changes setup: records returned in chunk are removed
      from setup.
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


@ray.remote
def resolve_transcripts_ray(primary_dat, secondary_dat, args):
    return resolve.resolve_transcripts(primary_dat, secondary_dat, args)


def collapse_list(dat_list, params):

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
                setup.append((dat_list[i], dat_list[i+1]))
            else:
                pass
        else:
           setup.append((dat_list[i], dat_list[i+1]))

    print(f"len(setup): {len(setup)}")

    ## resolve_transcripts(primary_dat, secondary_dat, params):
    if setup:
        futures = []
        for dat1, dat2 in setup:
            futures.append(resolve_transcripts_ray.remote(dat1, dat2, params))
        dat_collapse = ray.get(futures)
    else:
        dat_collapse = []

    if not last_paired:
        dat_collapse.append(dat_list[idx_last * 2])

    return dat_collapse


def update_xrefs(dat_list, xrefs):

    for dat in dat_list:
        for lost_id, kept_id in dat['filtered'].items():
            if lost_id in xrefs:
                print("\n" f"WARN: saw {lost_id} before\n")
            xrefs[lost_id] = kept_id


def process_gtfs(gtfs, params):

    setup = [(gtf, f"{i}.") for (gtf, i) in zip(gtfs, range(len(gtfs)))]
    dat = None
    xrefs = {}

    ## chunk: [(gtf1, 1), (gtf2, 2), ... ]; records removed from setup:
    chunk = pop_chunk(setup, params.nthreads)
    while chunk:

        futures = []
        for gtf, index in chunk:
            ## load_gtf(gtf_path, prefix, params); 
            ##   parses and resolves self:
            futures.append(load_gtf_ray.remote(gtf, f"{index}", params))

        dat_list = ray.get(futures)
        print(f"process_gtfs: load_gtf done at {output1.elapsed(params)} seconds")
        print(f"len(dat_list): {len(dat_list)}")
        update_xrefs(dat_list, xrefs)

        while len(dat_list) > 1:
            ## uses ray to parallelize calls to resolve_transcripts:
            dat_list = collapse_list(dat_list, params)
            print(
                f"process_gtfs: collapse_list done at "
                f"{output1.elapsed(params)} seconds"
            )
            print(f"len(dat_list): {len(dat_list)}")
            update_xrefs(dat_list, xrefs)

        if dat:
            ## resolve_transcripts(primary_dat, secondary_dat, args):
            dat = resolve.resolve_transcripts(dat, dat_list[0], params)    
            print(
                f"process_gtfs: resolve_transcripts done at "
                f"{output1.elapsed(params)} seconds"
            )
            update_xrefs([dat], xrefs)
        else:
            dat = dat_list[0]

        chunk = pop_chunk(setup, params.nthreads)

    dat['filtered'] = xrefs
    return dat


##########################################################################
## main:

if __name__ == '__main__':

    params = initial2.initialize()

    try:
        ray.init(num_cpus=params.nthreads, object_store_memory=params.memory)
        gtfs = read_gtf_list_file(params.gtf_list_file)
        dat = process_gtfs(gtfs, params)
    except Exception as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(31)

    print(f"main: process_gtfs done at {output1.elapsed(params)} seconds")

    output_filter = f"{params.output_prefix}.filter.tsv"
    try:
        with open(output_filter, 'w') as fh:
            fh.write("lost\tkept\n")
            for id1, id2 in dat['filtered'].items():
                fh.write(f"{id1}\t{id2}\n")
    except Exception as e:
        sys.stderr.write(f"ERROR: writing {output_filter}: {e}\n")
        sys.exit(33)

    print(
        f"main: finished writing {output_filter} at "
        f"{output1.elapsed(params)} seconds"
    )

    ids = assign_ids.assign_ids(dat, params)
    print(f"main: assign_ids done at {output1.elapsed(params)} seconds")

    ## [(src_gtf, old_id, new_id), ...]
    xref_table = output2.xrefs(dat['filtered'], ids, gtfs)
    print(f"main: xref_table done at {output1.elapsed(params)} seconds")

    output_xref = f"{params.output_prefix}.xref.tsv"
    try:
        with open(output_xref, 'w') as fh:
            header = '\t'.join([
                'source_gtf', 
                'old_id', 
                'transcript_id', 
                'gene_id'
            ])
            fh.write(f"{header}\n")
            for toks in xref_table:
                line = '\t'.join(toks)
                fh.write(f"{line}\n")
    except Exception as e:
        sys.stderr.write(f"ERROR: writing {output_xref}: {e}\n")
        sys.exit(33)

    print(
        f"main: finished writing {output_xref} at "
        f"{output1.elapsed(params)} seconds"
    )

    gtf_table = output2.dat2table(dat, ids)
    print(f"main: dat2table done at {output1.elapsed(params)} seconds")

    output_gtf = f"{params.output_prefix}.gtf"
    try:
        with open(output_gtf, 'w') as fh:
            for toks in gtf_table:
                toks = [str(tok) for tok in toks]
                line = '\t'.join(toks)
                fh.write(f"{line}\n")
    except Exception as e:
        sys.stderr.write(f"ERROR: writing {output_gtf}: {e}\n")
        sys.exit(35)

    print(
        f"main: finished writing {output_gtf} at "
        f"{output1.elapsed(params)} seconds"
    )

    sys.exit(0)


