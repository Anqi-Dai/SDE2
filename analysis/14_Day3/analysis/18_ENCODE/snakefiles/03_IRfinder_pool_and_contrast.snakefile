'''
This snakefile will do the following:
Use IFfinder to identify the differentially retend introns when comparing KD vs control
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

rule compare_KD_vs_ctrl:
    input:
        pooled_KD='{sample_dir}/IRfinder_result/Pooled_{gene}/IRFinder-IR-dir.txt',
        pooled_ctrl='{sample_dir}/IRfinder_result/Pooled_ctrl/IRFinder-IR-dir.txt',
        KD_1='{sample_dir}/IRfinder_result/{gene}-rep1/IRFinder-IR-dir.txt',
        KD_2='{sample_dir}/IRfinder_result/{gene}-rep2/IRFinder-IR-dir.txt',
        c11='{sample_dir}/IRfinder_result/ctrl-rep11/IRFinder-IR-dir.txt',
        c12='{sample_dir}/IRfinder_result/ctrl-rep12/IRFinder-IR-dir.txt',
        c2='{sample_dir}/IRfinder_result/ctrl-rep2/IRFinder-IR-dir.txt'
    output:
        '{sample_dir}/IRfinder_result/{gene}_vs_ctrl.txt'
    shell:
        '''
        /usr4/bs831/adai/bubhub-home/SDE2/analysis/05_IRfinder/IRfinder_software/IRFinder-1.2.5/bin/analysisWithLowReplicates.pl  \
            -A {input.pooled_KD} {input.KD_1} {input.KD_2}  \
            -B {input.pooled_ctrl} {input.c11} {input.c12}  {input.c2} > {output}
        '''

rule output:
    input:
        expand('{sample_dir}/IRfinder_result/{gene}_vs_ctrl.txt',
                sample_dir = sample_dir,
                gene = config['genes'])



