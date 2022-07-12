"""
assumes input gtf has both transcript and exon features and that both feature
  types have transcript_id, plus transcript features have gene_id;

this works w/ gencode and some refseq, but breaks on other refseq which do
  not have transcript features annotated; for those, would have to change
  the code in this module to only assume exon features w/ both transcript_id
  and gene_id; would then build transcript features from exon features; so
  necessary changes can be limited to this module.
"""

def parse_attributes(tok9, prefix=None):

    """
    Description: parse gtf2 attributes (9th column) for gene_id and transcript_id.

    Inputs:

      tok9: 9th column (attributes) from one row from input gtf file.

      prefix: for transcript_id; can be None.

    Output: dict w/ keys 'gene_id' and 'transcript_id'; set to None if missing,
    """

    toks = tok9.split(" ")
    toks = [tok_i.strip().rstrip(';').strip('"')  for tok_i in toks ]
 
    if prefix is None:
        prefix = ""

    if len(toks) < 4:
        raise Exception(f"Insufficient number of attribute tokens: {toks}")

    if len(toks) % 2:
        toks.pop()

    n = int(len(toks) / 2)
    attributes = { 'gene_id': None, 'transcript_id': None }

    for i in range(n):
        key = toks[2 * i]
        if key == 'gene_id':
            attributes['gene_id'] = f"{prefix}{toks[2 * i + 1]}"
        elif key == 'transcript_id':
            attributes['transcript_id'] = f"{prefix}{toks[2 * i + 1]}"

    return attributes


def register_transcript(toks, transcripts, locs=None, prefix=None):

    """
    Description: updates locs (if not None) and transcripts for one transcript
      record based on toks; transcripts defined by transcript_id; non-redundant
      collection saved.

    Inputs:
      toks is list of tokens, one per column from one row of input gtf file:

              0      1       2     3   4     5      6     7           8          9
        seqname source feature start end score strand frame [attributes] [comments]

        expected attributes e.g.: transcript_id "ENST00000456328.2";

      transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] ] }

      locs: { seq: [ [start, end, strand, transcript], ...], ..., }

    Output: primary transcript size; (1 + end - start);
    """

    attributes = parse_attributes(toks[8], prefix)

    if attributes['gene_id'] is None:
        raise Exception(f"register_transcript: could not find gene_id in {toks}.")

    if attributes['transcript_id'] is None:
        raise Exception(f"register_transcript: could not find transcript_id in {toks}.")

    transcripts[attributes['transcript_id']] = [
        attributes['gene_id'],       ## 0:gene_id
        toks[0],                     ## 1:seq_id
        toks[6],                     ## 2:strand
        toks[3],                     ## 3:start
        toks[4],                     ## 4:end
        []                           ## 5:exons
    ]                                ## python dumb

    if locs is not None:

        key = toks[0]                    ## seq

        val = [
            toks[3],                     ## 0:start
            toks[4],                     ## 1:end
            toks[1],                     ## 2:strand
            attributes['transcript_id']  ## 3:transcript_id
        ]

        if key not in locs:
            locs[key] = []

        locs[key].append(val)

    ## return length = 1 + end - start:
    return 1 + toks[4] - toks[3]


def register_exon(toks, exons, transcripts, prefix):

    """
    Description: Updates exons and transcripts for one exon record based on toks;
      exon defined by seq:strand:start:end; non-redundant collection saved.

    Inputs:
      toks is list of tokens, one per column from one row of input gtf file;

              0      1       2     3   4     5      6     7           8          9
        seqname source feature start end score strand frame [attributes] [comments]

        expected attributes e.g.: transcript_id "ENST00000456328.2";

      exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }

      transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exon] }

      prefix: for transcript_id and gene_id; can be None.

    Output: None
    """

    attributes = parse_attributes(toks[8], prefix)

    if attributes['transcript_id'] is None:
        raise Exception(f"register_exon: could not find transcript_id in {toks}.")

    ## seq:strand:start:end:
    key = f"{toks[0]}:{toks[6]}:{toks[3]}:{toks[4]}"

    exons[key] = [
        toks[0],             ## 0:seq
        toks[6],             ## 1:strand
        toks[3],             ## 2:start
        toks[4]              ## 3:end
    ]

    if attributes['transcript_id'] not in transcripts:
        raise Exception(
            f"register_exon: for exon '{key}', have not seen transcript_id "
            f"{attributes['transcript_id']} yet.")

    transcripts[attributes['transcript_id']][5].append(key)
 
    ## ensure exon_keys of transcript record stay sorted by (start, end),
    ##   so output order will be correct:

    def f_sort(exon_key):
        start = exons[exon_key][2]
        end = exons[exon_key][3]
        return (start, end)

    transcripts[attributes['transcript_id']][5].sort(key=f_sort)


def get_max_loc_length(loc):
    """
    """

    max_loc_length = 0
    return max_loc_length


def infer_transcripts(exons, prefix):

    """
    Inputs:
      exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }

    Output:
      transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] }

      locs: { seq: [ [start, end, strand, transcript], ...], ..., }
      transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] }
      exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }
      max_loc_length: length of longest encountered transcript

    """

    transcripts = {}
    return transcripts


def infer_locs(transcripts):
    """
    Inputs:
      transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] }

    Outputs:
      locs: { seq: [ [start, end, strand, transcript], ...], ..., }
    """

    locs = {}
    return locs


def parse_exons(fh, prefix):
    """
    """

    exons = {}
    return exons


def parse_gtf(fh, prefix):

    """
    Description: parses transcripts and exons from primary gtf;

    Inputs:

      fh: open readable filehandle to primary gtf;

      prefix: string to prepend to output transcript_ids;
        transcript_ids not changed if prefix is None.

    Outputs: dict w/ keys:
      locs: { seq: [ [start, end, strand, transcript], ...], ..., }
      transcripts: { 'transcript_id': [ gene_id, seq, strand, start, end, [exons] }
      exons: { 'seq:strand:start:end': [ seq, strand, start, end ] }
      max_tracnscript_length: length of longest encountered transcript
    """

    exons = {}
    transcripts = {}
    locs = {}
    max_loc_length = 0

    for line in fh:

        """
        toks is list of tokens, one per column from one row of input gtf file;
              0      1       2     3   4     5      6     7           8          9
        seqname source feature start end score strand frame [attributes] [comments]
        """

        toks = line.split('\t')
        toks = [ tok_i.strip() for tok_i in toks ]

        if len(toks) == 0:
            continue
        if toks[0].startswith('#'):
            continue
        if len(toks) < 9:
            raise Exception(f"parse_primary_gtf: fewer than 9 tokens on line: {line}.")

        ## start should be int:
        try:
            toks[3] = int(toks[3])
        except Exception as e:
            raise Exception(f"parse_primary_gtf: int() on '{toks[3]}': {e}")

        ## end should be int:
        try:
            toks[4] = int(toks[4])
        except Exception as e:
            raise Exception(f"parse_primary_gtf: int() on '{toks[4]}': {e}")

        if toks[2] == 'exon':
            register_exon(toks, exons, transcripts, prefix)
        else:
            continue

    max_loc_length = infer_transcripts(exons, transcripts, locs, prefix)

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

