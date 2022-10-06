# GTFMERGE
## Tools for merging GTF2.2 files.
#### Contact: mitch.kostich@jax.org
 
Installation:

```
## requires Python 3; developed and tested w/ Python 3.9.5

## install python dependencies (optionally put in virtual environment):
pip install -U argparse Cython ray setuptools

## download code to wherever you want it:
mkdir -p ~/opt
cd ~/opt
git clone https://github.com/TheJacksonLaboratory/gtfmerge.git

## optionally, build external modules for a modest speed gain:
cd gtfmerge
python setup.py build_ext --inplace

## in any case, put $HOME/opt/gtfmerge in your path or use full 
##   path to gtfmerge.py and mergegtfs.py commands. 
```

---

Currently there are two programs: `gtfmerge.py` and `mergegtfs.py`. Can use 
`gtfmerge.py` to merge a reference gtf with an isoform candidate gtf. Can use
`mergegtfs.py` to merge multiple isoform candidate gtfs. Neither program 
annotates gtfs, except `gtfmerge.py` will associate transcript identifiers
between the two input gtfs. Both programs will assign new 
gene_ids and transcript_ids, except `gtfmerge.py` will not if only given a
single input gtf. Clustering of transcripts into genes is 
controlled by cutoffs on the proportional overlap between exons and the 
proportion of all exons that overlap sufficiently. 

For eliminating redundancy in one GTF2.2 file or merging two GTF2.2 files,
usually use `gtfmerge.py`. This produces a cross-reference of deleted isoforms 
to kept isoforms. The `gtfmerge.py` program was primarily intended for merging 
a reference gtf with a consolidated isoform gtf:

```
$ ./gtfmerge.py -h
usage: gtfmerge.py [-h] [--tol_sj TOL_SJ] [--tol_tss TOL_TSS] [--tol_tts TOL_TTS] [--keep_all_primary] [--p_exon_overlap P_EXON_OVERLAP]
                   [--p_exons_overlap P_EXONS_OVERLAP] [--primary_prefix PRIMARY_PREFIX] [--secondary_prefix SECONDARY_PREFIX]
                   [--output_prefix OUTPUT_PREFIX] [--ids_from {1,2}]
                   primary_gtf [secondary_gtf]

Merges redundant transcripts from GTF2.2 formatted files. Merges primary_gtf with secondary_gtf, resulting in a non-redundant union gtf. Omitting
secondary_gtf results in non-redundant set from primary_gtf.

positional arguments:
  primary_gtf           GTF2.2 file containing primary isoform list
  secondary_gtf         GTF2.2 file with secondary isoform list (default: None)

optional arguments:
  -h, --help            show this help message and exit
  --tol_sj TOL_SJ       Tolerance (bp) for matching splice junction coordinates (default: 0)
  --tol_tss TOL_TSS     Tolerance (bp) for matching transcript start coordinates (default: 0)
  --tol_tts TOL_TTS     Tolerance (bp) for matching transcript end coordinates (default: 0)
  --keep_all_primary    Keep all primary transcripts (even redundant ones) (default: False)
  --p_exon_overlap P_EXON_OVERLAP
                        Minimum proportion overlap between two exons for gene matching (default: 0.5)
  --p_exons_overlap P_EXONS_OVERLAP
                        Minimum proportion of exons with overlaps needed for gene matching (default: 0.25)
  --primary_prefix PRIMARY_PREFIX
                        Prefix to be prepended to primary_gtf transcript_ids (default: None)
  --secondary_prefix SECONDARY_PREFIX
                        Prefix to be prepended to secondary_gtf transcript_ids (default: None)
  --output_prefix OUTPUT_PREFIX
                        Prefix for output files (default: union)
  --ids_from {1,2}      Source of transcript_ids and gene_ids for redundant transcripts; 1: primary_gtf; 2: secondary_gtf (default: 1)

Input GTF2.2 formatted files are required to have exon features with attributes that include 'gene_id' and 'transcript_id'. Output GTF2.2 formatted
file includes 'exon' and 'transcript' features, both with attributes (in order) 'gene_id', 'transcript_id', 'primary_id', 'secondary_id; where
primary_id is the transcript_id in primary_gtf, and secondary_id is the transcript_id in secondary_gtf; primary_id and/or secondary_id will be set
to empty string '' if the transcript is not present in the corresponding GTF file. Throws error if transcript_id in primary_gtf is also found in
secondary_gtf; set --primary_prefix and/or --secondary_prefix to fix this. When a transcript is found in both input files, but the exon boundaries
differ (within the specified tolerances), the output exon coordinates will be taken from primary_gtf. Also outputs a tab-delimited table cross-
referencing filtered transcript_ids to the kept transcript_id with which they were found to be redundant.
```

---

For merging multiple GTF2.2 files, use `mergegtfs.py`. It can handle many 
files at once. It was developed for merging isoform candidate lists from 
multiple samples. It assigns new gene_ids in a sensible and tunable way. Then 
transcript_ids are assigned sequentially based on gene_ids. One current 
limitation is that cross-reference table is not produced. This shortcoming 
will be addressed in a future update:

```
$ ./mergegtfs.py -h
usage: mergegtfs.py [-h] [--nthreads NTHREADS] [--memory_gb MEMORY_GB] [--tol_sj TOL_SJ] [--tol_tss TOL_TSS] [--tol_tts TOL_TTS] [--p_exon_overlap P_EXON_OVERLAP] [--p_exons_overlap P_EXONS_OVERLAP]
                    [--gene_prefix GENE_PREFIX] [--output_prefix OUTPUT_PREFIX]
                    gtf_list_file

Merges redundant transcripts from multiple GTF2.2 formatted files, resulting in a non-redundant union gtf.

positional arguments:
  gtf_list_file         File with list of GTF2.2 files (one path per line) to be merged

optional arguments:
  -h, --help            show this help message and exit
  --nthreads NTHREADS   Number of threads to use for parallel execution; even number more efficient (default: 1)
  --memory_gb MEMORY_GB
                        Amount of available RAM, in gigabytes (default: 8)
  --tol_sj TOL_SJ       Tolerance (bp) for matching splice junction coordinates (default: 0)
  --tol_tss TOL_TSS     Tolerance (bp) for matching transcript start coordinates (default: 0)
  --tol_tts TOL_TTS     Tolerance (bp) for matching transcript end coordinates (default: 0)
  --p_exon_overlap P_EXON_OVERLAP
                        Minimum proportion overlap between two exons needed for gene matching (default: 0.5)
  --p_exons_overlap P_EXONS_OVERLAP
                        Minimum proportion of overlapping exons needed for gene matching (default: 0.1)
  --gene_prefix GENE_PREFIX
                        Prefix for gene_ids and transcript_ids (default: PB.)
  --output_prefix OUTPUT_PREFIX
                        Prefix for output file names (default: union)

Input GTF2.2 formatted files are required to have exon features with attributes that include 'gene_id' and 'transcript_id'. Output GTF2.2 formatted file includes 'exon' and 'transcript' features, both
with attributes (in order) 'gene_id' and 'transcript_id'. New gene_ids will be assigned in accordance with --p_exon_overlap (min overlap for matching exons) and --p_exons_overlap (min proportion of
matched exons for matching gene_ids). New transcript_ids are assigned sequentially for each gene. New gene_ids and transcript_ids both have prefix specified by --gene_prefix.
```

