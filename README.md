'''
usage: gtfmerge.py [-h] [--tol_sj TOL_SJ] [--tol_tss TOL_TSS]
                   [--tol_tts TOL_TTS] [--ids_from {1,2}] [--keep_all_primary]
                   [--primary_prefix PRIMARY_PREFIX]
                   [--secondary_prefix SECONDARY_PREFIX]
                   primary_gtf [secondary_gtf]

Merges redundant transcripts from GTF2.2 formatted files. Merges primary_gtf
with secondary_gtf, resulting in a non-redundant union gtf. Omitting
secondary_gtf results in non-redundant set from primary_gtf.

positional arguments:
  primary_gtf           GTF2.2 file containing primary isoform list
  secondary_gtf         GTF2.2 file with secondary isoform list (default:
                        None)

optional arguments:
  -h, --help            show this help message and exit
  --tol_sj TOL_SJ       Tolerance (bp) for matching splice junction
                        coordinates (default: 0)
  --tol_tss TOL_TSS     Tolerance (bp) for matching transcript start
                        coordinates (default: 0)
  --tol_tts TOL_TTS     Tolerance (bp) for matching transcript end coordinates
                        (default: 0)
  --ids_from {1,2}      Source for transcript_ids and gene_ids; 1:
                        primary_gtf; 2: secondary_gtf (default: 1)
  --keep_all_primary    Keep all primary transcripts (even redundant ones)
                        (default: False)
  --primary_prefix PRIMARY_PREFIX
                        Prefix to be prepended to primary_gtf transcript_ids
                        (default: None)
  --secondary_prefix SECONDARY_PREFIX
                        Prefix to be prepended to secondary_gtf transcript_ids
                        (default: None)

Input GTF2.2 formatted files are required to have exon features with
attributes that include 'gene_id' and 'transcript_id'. Output GTF2.2 formatted
file includes 'exon' and 'transcript' features, both with attributes (in
order) 'gene_id', 'transcript_id', 'primary_id', 'secondary_id; where
primary_id is the transcript_id in primary_gtf, and secondary_id is the
transcript_id in secondary_gtf; primary_id and/or secondary_id will be set to
empty string '' if the transcript is not present in the corresponding GTF
file. Error if transcript_id in primary_gtf is also found in secondary_gtf;
set --primary_prefix and/or --secondary_prefix to fix this. When a transcript
is found in both input files, but the exon boundaries differ (within the
specified tolerances), the output exon coordinates will be taken from
primary_gtf.
'''
