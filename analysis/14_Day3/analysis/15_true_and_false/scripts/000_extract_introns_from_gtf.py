#!/usr/bin/env python
'''
This script accepts a (not necessarily sorted) gtf file as its only command line argument,
extracts and groups all 'exon' features by transcript_id,
computes distinct intronic coordinates on the gene_id level,
and writes out the results sorted by chromosome and start coordinate
to stdout in bed format, e.g.:

chr1	12057	12179	+	ENSG00000223972.5
chr1	12227	12613	+	ENSG00000223972.5
'''
from collections import defaultdict
from csv import reader, writer
import re
import sys
fn = sys.argv[1]

transcripts = defaultdict(list)

tid_regex = re.compile('transcript_id "([^"]+)"')
with open(fn) as f : 
    f = reader(f,delimiter='\t')

    for rec in f : 
        if rec[0].startswith('#') : 
            continue
      
        if rec[2] == 'exon' : 
            tid = tid_regex.search(rec[-1]).group(1)
            transcripts[tid].append(rec)

genes = defaultdict(list)
gid_regex = re.compile('gene_id "([^"]+)"')
for tid, exons in transcripts.items() : 
    gid = gid_regex.search(exons[0][-1]).group(1)
    exons = sorted(exons,key=lambda x: int(x[3]))
    for e1, e2 in zip(exons,exons[1:]) : 
        genes[gid].append((exons[0][0],int(e1[4]),int(e2[3]),exons[0][6],gid))

all_introns = []
for gid, introns in genes.items() : 
    all_introns.extend(list(set(introns)))

out_f = writer(sys.stdout,delimiter='\t')
out_f.writerows(sorted(all_introns))
