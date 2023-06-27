from pipeline_util import *

localrules: all_reads_coverage_results

######## Step: Visualization: Reads coverage ##################
rule all_reads_coverage_results:
    input:
         bws=expand('bw/{sample}.bw', sample=aligned_names(config, FASTQ_PATHS, BAMS_PATHS)),

rule index_bams:
    input: '{anywhere}/{sample}.bam'
    output: '{anywhere}/{sample}.bam.bai'
    log: 'logs/bams_index/{anywhere}/{sample}.bam.bai.log'
    conda: '../envs/bio.env.yaml'
    shell: 'samtools index {input} {output} &>> {log}'

rule bam2bw:
    input:
         bam=f"{config['bams_dir']}/{{filename}}.bam",
         bai=f"{config['bams_dir']}/{{filename}}.bam.bai"
    output: 'bw/{filename, [^/]*}.bw'
    log: 'logs/bw/{filename}.log'
    conda: '../envs/deeptools.env.yaml'
    threads: 4
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell: 'bamCoverage -b {input.bam} -p {threads} -o {output} &>> {log}'
