'''
This snakefile will do the following:
* pooling all KD as experiment and run the IRfinder
'''

from glob import glob
import os

# -----------------------------------------------------------

sample_dir = '../../samples/replicate1'

sample_names = json.load(open('../samples.json')).keys()

DNA=os.path.realpath('ref/v27/genome.fa')
GTF=os.path.realpath('ref/v27/transcripts.gtf')

IRfinder_path = 'IRfinder_software/IRFinder-1.2.5/bin/IRFinder'

# -----------------------------------------------------------

rule all:
    input:
        sample_dir + '/IRfinder_result/RFKD_vs_ctrl.txt'

def get_RF14_3replicates_unsorted_bam(wildcards):
    return glob('{}/IRfinder_result/RF-14*/Unsorted.bam'.format(wildcards.sample_dir))

def get_RF89_3replicates_unsorted_bam(wildcards):
    return glob('{}/IRfinder_result/RF-89*/Unsorted.bam'.format(wildcards.sample_dir))

def get_control_3replicates_unsorted_bam(wildcards):
    return glob('{}/IRfinder_result/RF-Control*/Unsorted.bam'.format(wildcards.sample_dir))

rule concate_all_experiment:
    input:
        RF14=get_RF14_3replicates_unsorted_bam,
        RF89=get_RF89_3replicates_unsorted_bam
    output:
        outPooled='{sample_dir}/IRfinder_result/concat_all_KD_unsorted.bam'
    threads:
        4
    shell:
        '''
        samtools cat {input.RF14} {input.RF89} -o {output.outPooled};
        '''

rule pooled_KD:
    input:
        bam='{sample_dir}/IRfinder_result/concat_all_KD_unsorted.bam',
        ref='REF/Human-hg38',
        IRpath=IRfinder_path 
    output:
        '{sample_dir}/IRfinder_result/Pooled_KD/IRFinder-IR-dir.txt'
    params:
        prefix='{sample_dir}/IRfinder_result/Pooled_KD/'
    shell:
        '''
        export LD_LIBRARY_PATH=/projectnb/bubhub/conda_root/user_conda/adai/envs/flynn_altstatus/lib;
        {input.IRpath} -m BAM -r {input.ref} -d  {params.prefix} {input.bam} > {output}
        '''    



rule compare_IR_KDvsctrl:
    input:
        pooled_KD='{sample_dir}/IRfinder_result/Pooled_KD/IRFinder-IR-dir.txt',
        pooled_ctrl='{sample_dir}/IRfinder_result/Pooled_ctrl/IRFinder-IR-dir.txt',
        KD_1='{sample_dir}/IRfinder_result/RF-14-1/IRFinder-IR-dir.txt',
        KD_2='{sample_dir}/IRfinder_result/RF-14-2/IRFinder-IR-dir.txt',
        KD_3='{sample_dir}/IRfinder_result/RF-14-3/IRFinder-IR-dir.txt',
        KD_4='{sample_dir}/IRfinder_result/RF-89-1/IRFinder-IR-dir.txt',
        KD_5='{sample_dir}/IRfinder_result/RF-89-2/IRFinder-IR-dir.txt',
        KD_6='{sample_dir}/IRfinder_result/RF-89-3/IRFinder-IR-dir.txt',
        c1='{sample_dir}/IRfinder_result/RF-Control-1/IRFinder-IR-dir.txt',
        c2='{sample_dir}/IRfinder_result/RF-Control-2/IRFinder-IR-dir.txt',
        c3='{sample_dir}/IRfinder_result/RF-Control-3/IRFinder-IR-dir.txt'
    output:
        '{sample_dir}/IRfinder_result/RFKD_vs_ctrl.txt'
    shell:
        '''
        IRfinder_software/IRFinder-1.2.5/bin/analysisWithLowReplicates.pl  \
            -A {input.pooled_KD} {input.KD_1} {input.KD_2} {input.KD_3} {input.KD_4} {input.KD_5}  {input.KD_6}\
            -B {input.pooled_ctrl} {input.c1} {input.c2}  {input.c3}  > {output}
        '''


