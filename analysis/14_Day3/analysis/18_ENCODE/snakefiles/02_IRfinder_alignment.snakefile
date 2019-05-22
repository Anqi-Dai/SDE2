'''
This snakefile will do the following:

* Map the reads to the reference and produce the bam file
* Concat the unsorted bam file
'''

from glob import glob
import os

# -----------------------------------------------------------

sample_dir = '/usr4/bs831/adai/adai/ENCODE_KD_samples'

sample_names = json.load(open('../samples.json')).keys()

IRfinder_path = '/usr4/bs831/adai/bubhub-home/SDE2/analysis/05_IRfinder/IRfinder_software/IRFinder-1.2.5/bin/IRFinder'
Ref_dir = '/usr4/bs831/adai/bubhub-home/SDE2/analysis/05_IRfinder/REF/Human-hg38'

configfile : 'IRfinder_config.yml' 
# -----------------------------------------------------------

rule all:
    input:
        expand('{sample_dir}/IRfinder_result/{sample_names}/Unsorted.bam',sample_dir = sample_dir, sample_names = sample_names)
   

# Align the reads to the reference through IRfinder, cuz I need the IRfinder part of the alignment result for the later IR differential analysis.
#   so the file structure has to be like specified in the documentation

rule Align_reads:
    input:
        R1='{sample_dir}/{sample_names}_R1_paired.fastq.gz',
        R2='{sample_dir}/{sample_names}_R2_paired.fastq.gz',
        ref=Ref_dir,
        IRpath=IRfinder_path
    output:
        bam='{sample_dir}/IRfinder_result/{sample_names}/Unsorted.bam',
        txt='{sample_dir}/IRfinder_result/{sample_names}/IRFinder-IR-dir.txt'
    params:
        prefix='{sample_dir}/IRfinder_result/{sample_names}/'
    threads:
        16
    shell:
        '''
        export LD_LIBRARY_PATH=/projectnb/bubhub/conda_root/user_conda/adai/envs/flynn_altstatus/lib;
        {input.IRpath} -d  {params.prefix}  -s LoadAndKeep \
         -u -t {threads}  -r {input.ref} {input.R1} {input.R2}
        '''


# Concat the replicates together
# wow I can use the config file to actually do a for loop which is making the code too much easier to write and cleaner.

rule output:
        input:
            expand('{sample_dir}/IRfinder_result/{gene}_concat_rep_unsorted.bam', sample_dir = sample_dir,  gene = config['genes'])

def gather_KD_bio_replicates(wildcards):
    return expand('{sample_dir}/IRfinder_result/{gene}-rep{replicated_num}/Unsorted.bam',
                    sample_dir= wildcards.sample_dir,
                    gene = wildcards.gene,
                    replicated_num = [1,2])

rule concat_replicates:
    input:
        gather_KD_bio_replicates
    output:
        '{sample_dir}/IRfinder_result/{gene}_concat_rep_unsorted.bam'  
    threads:
        4
    shell:
        '''
        samtools cat  {input} -o {output}
        '''

# concat the control samples together

def gather_control_replicates(wildcards):
    return expand('{sample_dir}/IRfinder_result/ctrl-rep{replicated_num}/Unsorted.bam',
                    sample_dir= wildcards.sample_dir,
                    replicated_num = [11,12,2 ])


rule concat_control:
    input:
        gather_control_replicates
    output:
        '{sample_dir}/IRfinder_result/ctrl_concat_rep_unsorted.bam'
    threads:
        4
    shell:
        '''
        samtools cat  {input} -o {output}
        '''

rule output2:
    input: sample_dir + '/IRfinder_result/ctrl_concat_rep_unsorted.bam'
