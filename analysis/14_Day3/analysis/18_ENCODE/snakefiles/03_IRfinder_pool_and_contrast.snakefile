'''
This snakefile will do the following:
Use IFfinder to identify the differentially retend introns when comparing KD vs control

* Build reference for IRfinder
* Map the reads to the reference and produce the bam file
* Concat the unsorted bam file
* Analyze using analysisWithLowReplicates.pl  
'''

from glob import glob
import os

# -----------------------------------------------------------
configfile: 'IRfinder_config.yml'

sample_dir = '/usr4/bs831/adai/adai/ENCODE_KD_samples'

sample_names = json.load(open('../samples.json')).keys()

IRfinder_path = '/usr4/bs831/adai/bubhub-home/SDE2/analysis/05_IRfinder/IRfinder_software/IRFinder-1.2.5/bin/IRFinder'

Ref_dir = '/usr4/bs831/adai/bubhub-home/SDE2/analysis/05_IRfinder/REF/Human-hg38'

# -----------------------------------------------------------

rule output:
    input:
        expand('{sample_dir}/IRfinder_result/Pooled_{each}/IRFinder-IR-dir.txt',
                sample_dir = sample_dir,
                each = config['all'])


rule pooling_IR:
    input:
        bam='{sample_dir}/IRfinder_result/{each}_concat_rep_unsorted.bam',
        ref=Ref_dir,
        IRpath=IRfinder_path 
    output:
        '{sample_dir}/IRfinder_result/Pooled_{each}/IRFinder-IR-dir.txt'
    threads:
        12
    params:
        prefix='{sample_dir}/IRfinder_result/Pooled_{each}/'
    shell:
        '''
        export LD_LIBRARY_PATH=/projectnb/bubhub/conda_root/user_conda/adai/envs/flynn_altstatus/lib;
        {input.IRpath} -m BAM -r {input.ref} -d  {params.prefix} {input.bam} > {output}
        '''    

# Calling the comparison script to get the differential IR results

rule compare_IR_14vsctrl:
    input:
        pooled_14='{sample_dir}/IRfinder_result/Pooled_14/IRFinder-IR-dir.txt',
        pooled_ctrl='{sample_dir}/IRfinder_result/Pooled_ctrl/IRFinder-IR-dir.txt',
        RF14_1='{sample_dir}/IRfinder_result/RF-3-14-1/IRFinder-IR-dir.txt',
        RF14_2='{sample_dir}/IRfinder_result/RF-3-14-2/IRFinder-IR-dir.txt',
        RF14_3='{sample_dir}/IRfinder_result/RF-3-14-3/IRFinder-IR-dir.txt',
        c1='{sample_dir}/IRfinder_result/RF-3-M-1/IRFinder-IR-dir.txt',
        c2='{sample_dir}/IRfinder_result/RF-3-M-2/IRFinder-IR-dir.txt',
        c3='{sample_dir}/IRfinder_result/RF-3-M-3/IRFinder-IR-dir.txt'
    output:
        '{sample_dir}/IRfinder_result/RF14_vs_ctrl.txt'
    shell:
        '''
        ../../../05_IRfinder/IRfinder_software/IRFinder-1.2.5/bin/analysisWithLowReplicates.pl  \
            -A {input.pooled_14} {input.RF14_1} {input.RF14_2} {input.RF14_3}  \
            -B {input.pooled_ctrl} {input.c1} {input.c2}  {input.c3} > {output}
        '''

