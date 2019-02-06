'''
This Snakefile will generate detk report on the raw counts of SDE2 data(the raw counts from csvgather)
for educational purpose of the de_toolkit workflow.
'''

rule normalize_raw:
    input:
        'data/raw_salmon_cts_replicate1_matrix.csv'
    output:
        'output/norm_salmon_cts_replicate1_matrix_from_raw.csv'
    shell:
        'detk-norm deseq2  {input} -o {output}'

rule stats:
    input:
        counts = 'data/raw_salmon_cts_replicate1_matrix.csv',
        pheno= 'data/pheno_data.csv'
    output:
        'output/raw_summary_stats.csv'
    shell:
        'detk-stats summary --log {input.counts} -o {output} --column-data={input.pheno}'
    
rule generate:
    input:
        'output/raw_summary_stats.csv'
    output:
        'detk_report/detk_report.html'
    shell:
        'detk-report generate --dev'


rule test_stats:
    input:
        counts = 'output/norm_salmon_cts_replicate1_matrix_from_raw.csv',
        pheno= 'data/pheno_data.csv'
    output:
        'test/norm_salmon_stats_replicate1_matrix_from_raw.csv'
    shell:
        'detk-stats summary --log {input.counts} -o {output} --column-data={input.pheno}'


rule test_generate:
    input:
        'test/norm_salmon_stats_replicate1_matrix_from_raw.csv'
    output:
        'detk_report_test/detk_report.html'
    shell:
        'detk-report generate --dev'
