from pipeline_util import *

localrules: all_reads_bw_results

######## Step: Visualization: Reads baws ##################
rule all_reads_bw_results:
    input:
         bws=expand('bw/{sample}.bw', sample=aligned_names(config, FASTQ_PATHS, BAMS_PATHS)) \
             if bool(config['bw']) else []

rule all_tags_bw_results:
    input:
        bws=expand('tagsbw/{sample}.bw', sample=aligned_names(config, FASTQ_PATHS, BAMS_PATHS)) \
            if bool(config['tagsbw']) else []


rule index_bams:
    input: '{anywhere}/{sample}.bam'
    output: '{anywhere}/{sample}.bam.bai'
    log: 'logs/bams_index/{anywhere}/{sample}.bam.bai.log'
    conda: '../envs/bio.env.yaml'
    shell: 'samtools index {input} {output} &> {log}'

rule bam2bw:
    input:
         bam=f"{config['bams_dir']}/{{sample}}.bam",
         bai=f"{config['bams_dir']}/{{sample}}.bam.bai"
    output: 'bw/{sample, [^/]*}.bw'
    log: 'logs/bw/{sample}.log'
    conda: '../envs/deeptools.env.yaml'
    threads: 4
    params:
        bamCoverage_params=config['bamCoverage_params'],
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell: 'bamCoverage -b {input.bam} -p {threads} -o {output} {params.bamCoverage_params} &> {log}'

rule bam2tagsbw:
    input:
        deduplicated_bam="deduplicated/{sample}.bam",
        chrom_sizes=rules.download_chrom_sizes.output
    output: 'tagsbw/{sample, [^/]*}.bw'
    log: 'logs/tagsbw/{sample}.log'
    conda: '../envs/bio.env.yaml'
    threads: 1
    params:
        fragment=config['tags2bw_fragment'],
    resources:
        mem = 6, mem_ram = 4,
        time = 60 * 120
    shell: 'bash {workflow.basedir}/scripts/bam2tagsbw.sh \
            {input.deduplicated_bam} {params.fragment} {input.chrom_sizes} {output} &> {log}'

rule bam_sorted:
    input: f"{config['bams_dir']}/{{sample}}.bam"
    output: temp('sorted/{sample}.bam')
    log: 'logs/bam_sorted/{sample}.log'
    conda: '../envs/bio.env.yaml'
    threads: 1
    resources:
        mem=6, mem_ram=4,
        time=60 * 120
    shell:
        'samtools sort {input} -o {output}  &> {log}'

rule bam_deduplicated:
    input: "sorted/{sample}.bam"
    output:
        bam=temp("deduplicated/{sample}.bam"),
        metrics=temp("qc/picard/{sample}.txt")
    log: "logs/bam_deduplicated/{sample}.log"
    threads: 1
    resources:
        mem = 6, mem_ram = 4,
        time = 60 * 120
    params: "REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=SILENT"
    wrapper: "0.31.1/bio/picard/markduplicates"
