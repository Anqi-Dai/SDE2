'''
# This Snakefile is the pipeline for the processes:
Trimming the fastq files of the ENCODE data (8 gene kd experiments)
'''

import glob
import json
import os

# Globals ---------------------------------------------------------------------

# Sample names
FILES = json.load(open('../samples.json'))
sample_names = sorted(FILES.keys())

# The samples directory
sample_dir = '/usr4/bs831/adai/adai/ENCODE_KD_samples'


# Rules -----------------------------------------------------------------------

rule all:
    input:
        expand('{sample_dir}/{sample_names}_R1_paired.fastq.gz',sample_dir = sample_dir, sample_names = sample_names)

#trimmomatic to trim the reads

rule trimmomatic:
    input:
        R1 = '{sample_dir}/{partname}_R1.fastq.gz',
        R2 = '{sample_dir}/{partname}_R2.fastq.gz'
    output:
        paired1 = '{sample_dir}/{partname}_R1_paired.fastq.gz',
        unpaired1 = '{sample_dir}/{partname}_R1_unpaired.fastq.gz',
        paired2 = '{sample_dir}/{partname}_R2_paired.fastq.gz',
        unpaired2 = '{sample_dir}/{partname}_R2_unpaired.fastq.gz',
        res_txt = '{sample_dir}/{partname}.txt'
    threads:
        8
    shell:
        """ java -jar /projectnb/bubhub/conda_root/user_conda/adai/envs/flynn_altstatus/share/trimmomatic-0.38-1/trimmomatic.jar \
            PE {input.R1} {input.R2} {output.paired1} {output.unpaired1} {output.paired2} {output.unpaired2} \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
            -threads {threads} \
            2> {output.res_txt}; """

