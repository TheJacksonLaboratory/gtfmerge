#!/usr/bin/env python3

## system:
import ray
import sys
import time

## local:
import assign_ids
import initial2
import gtf
import output
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


def process_gtfs(params):

    gtfs = read_gtf_list_file(params.gtf_list_file)
    setup = [(gtf, f"{i}.") for (gtf, i) in zip(gtfs, range(len(gtfs)))]
    dat = None

    chunk = pop_chunk(setup, params.nthreads)
    while chunk:

        futures = []
        for gtf, index in chunk:
            ## load_gtf(gtf_path, prefix, params); 
            ##   parses and resolves self:
            futures.append(load_gtf_ray.remote(gtf, f"{index}", params))

        dat_list = ray.get(futures)
        print(f"process_gtfs: load_gtf done at {output.elapsed(params)} seconds")
        print(f"len(dat_list): {len(dat_list)}")

        while len(dat_list) > 1:
            ## uses ray to parallelize calls to resolve_transcripts:
            dat_list = collapse_list(dat_list, params)
            print(
                f"process_gtfs: collapse_list done at "
                f"{output.elapsed(params)} seconds"
            )
            print(f"len(dat_list): {len(dat_list)}")

        if dat:
            ## resolve_transcripts(primary_dat, secondary_dat, args):
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
        ray.init(num_cpus=params.nthreads, object_store_memory=params.memory)
        dat = process_gtfs(params)
    except Exception as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(31)

    print(f"main: process_gtfs done at {output.elapsed(params)} seconds")

    ids = assign_ids.resolve_ids(dat, params)
    print(f"main: resolve_ids done at {output.elapsed(params)} seconds")

    gtf_table = output2.dat2table(dat, ids)
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


