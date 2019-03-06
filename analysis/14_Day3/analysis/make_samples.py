#!/usr/bin/env python3
'''
make a samples.json file with sample names and file names
of the SDE2 samples.
'''

import json
from glob import glob

fastqs = glob('/restricted/projectnb/flynngrp/adai/Day3_SDE2/*_001.fastq.gz')
FILES = {}

# Change this line to extract a sample name from each filename.
SAMPLES = [fastq.split('/')[-1].split('_')[0] for fastq in fastqs]
# use set to retrieve unique sample names
SAMPLES = set(SAMPLES)

for sample in SAMPLES:
    # Change 'R1' and 'R2' to match the way your mate pairs are marked.
        FILES[sample] = {}
        FILES[sample]['R1'] = sorted([fastq for fastq in fastqs if sample in fastq and 'R1' in fastq])
        FILES[sample]['R2'] = sorted([fastq for fastq in fastqs if sample in fastq and 'R2' in fastq])

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples.json', 'w').writelines(js)
