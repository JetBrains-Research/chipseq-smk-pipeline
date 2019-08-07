from pipeline_util import *

localrules: step5_reads_coverage_results

######## Step: Visualization: Reads coverage ##################
rule step5_reads_coverage_results:
    input:
         bws=expand('bw/{sample}.bw', sample=fastq_aligned_names(config)),

rule index_bams:
    input: '{anywhere}/{sample}.bam'
    output: '{anywhere}/{sample}.bam.bai'
    log: 'logs/bams_index/{anywhere}/{sample}.bam.bai.log'
    wrapper: '0.36.0/bio/samtools/index'

rule bam2bw:
    input:
         bam='bams/{filename}.bam',
         bai='bams/{filename}.bam.bai'
    output: 'bw/{filename, [^/]*}.bw'
    log: 'logs/bw/{filename}.log'

    conda: '../envs/deeptools.env.yaml'
    threads: 4
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell: 'bamCoverage -b {input.bam} -p {threads} -o {output}'
